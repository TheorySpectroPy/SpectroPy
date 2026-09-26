from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from string import ascii_lowercase
import os
import numpy as np
import yaml

from process_symmetry import read_site_symmetries, read_symmetry_file
from spectropy_config import read_settings
from spectropy_structure import Structure, read_structure, write_structure


@dataclass(frozen=True)
class Displacement:
    atom_index: int
    fractional_vector: np.ndarray
    suffix: str

    @property
    def name(self) -> str:
        return f"atom{self.atom_index}{self.suffix}"


def atom_labels(symbols: list[str]) -> list[str]:
    counts: dict[str, int] = defaultdict(int)
    labels: list[str] = []
    for symbol in symbols:
        counts[symbol] += 1
        labels.append(f"{symbol}{counts[symbol]}")
    return labels


def read_displacements(path: str | Path = "displacements.dat") -> list[Displacement]:
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Incomplete displacement file: {path}")
    count = int(lines[1].split()[0])
    displacements: list[Displacement] = []
    counters: dict[int, int] = defaultdict(int)
    for line in lines[2:2 + count]:
        values = line.split()
        if len(values) < 5:
            raise ValueError(f"Invalid displacement record: {line}")
        atom_index = int(values[1])
        suffix_index = counters[atom_index]
        if suffix_index >= len(ascii_lowercase):
            raise ValueError(f"Too many displacements for atom {atom_index}")
        displacements.append(Displacement(atom_index, np.array(values[2:5], dtype=float), ascii_lowercase[suffix_index]))
        counters[atom_index] += 1
    if len(displacements) != count:
        raise ValueError(f"Expected {count} displacements in {path}")
    return displacements


def read_displacement_mode(path: str | Path = "displacements.dat") -> str:
    header = Path(path).read_text().splitlines()[0].lower()
    if "fractional" in header:
        return "fractional"
    if "site-symmetry-minimal" in header:
        return "minimal"
    if "symmetry-inequivalent" in header:
        return "atoms"
    if header == "displacements for vasp. system":
        return "full"
    raise ValueError(f"Cannot determine displacement mode from {path}")


def write_displacements(
    path: str | Path,
    displacements: list[Displacement],
    structure: Structure,
    header: str = "displacements for VASP. System",
) -> None:
    labels = atom_labels(structure.atom_symbols)
    unique_atoms = list(dict.fromkeys(displacement.atom_index for displacement in displacements))
    with open(path, "w") as output:
        output.write(f"{header}\n")
        output.write(f"{len(displacements):6d}   Number of displacements\n")
        for displacement in displacements:
            vector = displacement.fractional_vector
            label = labels[displacement.atom_index - 1]
            output.write(
                f"{label:<5s}{displacement.atom_index:6d}"
                f"{vector[0]:18.12f}{vector[1]:18.12f}{vector[2]:18.12f}\n"
            )
        output.write(f"{len(unique_atoms):6d}{structure.natoms:6d}     Number of atoms in SC\n")
        for atom_index in unique_atoms:
            label = labels[atom_index - 1]
            output.write(f"{label:<5s}{atom_index:6d}\n")


def _amplitude(input_path: str | Path) -> np.ndarray:
    """Per-axis Cartesian amplitudes; one input value applies to all three axes."""
    value = read_settings(input_path).displacement_amplitude
    return np.array(value if value is not None else (0.03, 0.03, 0.03))


def _cartesian_vectors(amplitude: np.ndarray) -> np.ndarray:
    vectors = np.zeros((6, 3))
    for axis in range(3):
        vectors[2 * axis, axis] = amplitude[axis]
        vectors[2 * axis + 1, axis] = -amplitude[axis]
    return vectors


def _representatives(path: str | Path) -> list[int]:
    data = yaml.safe_load(Path(path).read_text()) or {}
    mapping = data.get("atom_mapping")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{path} does not contain Phonopy atom_mapping data")
    representatives = sorted({int(value) for value in mapping.values()})
    return representatives


# ---------------------------------------------------------------------------
# Cartesian minimal set ("minimal" mode)
# ---------------------------------------------------------------------------

def _frac_to_cart_rotation(rotation: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """Cartesian-basis form of a fractional-basis rotation.

    Matches reconstruct_dielectric_derivatives.frac_to_cart_rotation so the
    generator and the derivative reconstruction share one convention.
    """
    return lattice.T @ rotation @ np.linalg.inv(lattice.T)


def _get_displacement(rotations, directions, tol: float = 1e-8):
    """Minimal subset of ``directions`` whose site-symmetry orbits span 3D.

    Tries one direction, then two, then all three (Phonopy's one/two/three
    reduction). ``rotations`` are Cartesian site-symmetry rotations, so the
    orbit of a direction is ``{R @ d}``.
    """
    for size in (1, 2, 3):
        for subset in combinations(directions, size):
            orbit = [rotation @ direction for direction in subset for rotation in rotations]
            if np.linalg.matrix_rank(np.array(orbit), tol=tol) == 3:
                return list(subset)
    return list(directions)


def _need_minus(direction, rotations, tol: float = 1e-8) -> bool:
    """Add -d unless -d already lies in the site-symmetry orbit of d.

    Mirrors Phonopy's is_minus_displacement. When the negative is not needed,
    reconstruct_dielectric_derivatives.synthesize_missing_sign recovers it from
    the same site symmetry.
    """
    return not any(
        np.allclose(rotation @ direction, -direction, atol=tol) for rotation in rotations
    )


def _minimal_entries(structure: Structure, symmetry_path: str | Path, amplitude: np.ndarray, is_plusminus="auto") -> list[tuple[int, np.ndarray]]:
    """Cartesian symmetry-minimal displacement set.

    Uses Phonopy's one/two/three reduction and +/- rule, but with the Cartesian
    axes {x, y, z} as candidates, so every emitted direction is axis-aligned.
    The rank-3 dielectric derivative transforms covariantly under the site
    symmetry, so probing along a Cartesian axis is always possible and never
    needs more displacements than Phonopy's lattice-basis set. Site symmetries
    are read from the same ``symmetry`` file the reconstruction uses.
    """
    lattice = structure.lattice
    _, _, atom_mapping = read_symmetry_file(symmetry_path)
    site_symmetries = read_site_symmetries(symmetry_path)
    representatives = sorted(set(atom_mapping.values()) | set(site_symmetries))
    directions = np.eye(3)  # Cartesian x, y, z
    inverse_lattice = np.linalg.inv(lattice)

    vectors: list[tuple[int, np.ndarray]] = []
    for atom in representatives:
        fallback = np.eye(3)[None, :, :]
        fractional_rotations = site_symmetries.get(atom, {}).get("rotations", fallback)
        rotations = [_frac_to_cart_rotation(rotation, lattice) for rotation in fractional_rotations]
        for direction in _get_displacement(rotations, directions):
            vector = direction * amplitude
            vectors.append((atom, vector))
            if is_plusminus is True or (is_plusminus == "auto" and _need_minus(direction, rotations)):
                vectors.append((atom, -vector))
    return [(atom, vector @ inverse_lattice) for atom, vector in vectors]


# ---------------------------------------------------------------------------
# Fractional-basis minimal set ("fractional" mode) - Phonopy's exact algorithm
# ---------------------------------------------------------------------------

# Candidate displacement directions in the lattice (fractional) basis, in
# Phonopy's ``directions_diag`` order.  The one/two/three reduction below tries
# them in this exact order, so the emitted set matches Phonopy's.
DIRECTIONS_DIAG = np.array(
    [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
        [1, -1, 0],
        [1, 0, -1],
        [0, 1, -1],
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [-1, 1, 1],
    ],
    dtype=float,
)


def _det(a, b, c):
    """Determinant of the 3x3 matrix whose rows are a, b, c (Phonopy form)."""
    return (
        a[0] * b[1] * c[2]
        - a[0] * b[2] * c[1]
        + a[1] * b[2] * c[0]
        - a[1] * b[0] * c[2]
        + a[2] * b[0] * c[1]
        - a[2] * b[1] * c[0]
    )


def _rotated(direction, rotations):
    """All site-symmetry images of a fractional direction (Phonopy convention)."""
    return [direction @ rotation.T for rotation in rotations]


def _get_displacement_one(rotations, directions):
    """One lattice direction whose site-symmetry orbit spans 3D."""
    for direction in directions:
        images = _rotated(direction, rotations)
        for i in range(len(rotations)):
            for j in range(i + 1, len(rotations)):
                if _det(direction, images[i], images[j]) != 0:
                    return [direction]
    return None


def _get_displacement_two(rotations, directions):
    """Two lattice directions whose site-symmetry orbit spans 3D."""
    for direction in directions:
        images = _rotated(direction, rotations)
        for i in range(len(rotations)):
            for second in directions:
                if _det(direction, images[i], second) != 0:
                    return [direction, second]
    return None


def _get_fractional_displacement(rotations, directions):
    """Minimal set of lattice directions: one, then two, then three (Phonopy)."""
    result = _get_displacement_one(rotations, directions)
    if result is not None:
        return result
    result = _get_displacement_two(rotations, directions)
    if result is not None:
        return result
    return [directions[0], directions[1], directions[2]]


def _need_minus_fractional(direction, rotations) -> bool:
    """Add -d unless a site-symmetry operation maps d onto -d (Phonopy)."""
    for rotation in rotations:
        if not (direction @ rotation.T + direction).any():
            return False
    return True


def _fractional_entries(structure: Structure, symmetry_path: str | Path, amplitude: np.ndarray, is_plusminus="auto") -> list[tuple[int, np.ndarray]]:
    """Phonopy-equivalent fractional-basis site-symmetry-minimal displacement set.

    Replicates ``phonopy.harmonic.displacement.get_least_displacements`` with
    ``is_diagonal=True``: for every symmetry-inequivalent atom (an atom that
    maps onto itself in ``atom_mapping``) the site symmetry expands one, two,
    or three lattice-vector directions into a set spanning 3D, and each
    direction is paired with its negative unless the site symmetry already
    maps the direction onto its negative.  The selected directions are then
    converted to Cartesian and scaled to a uniform magnitude (Phonopy's single
    ``distance``; if the per-axis amplitudes differ, the first value is used).
    """
    lattice = structure.lattice
    _, _, atom_mapping = read_symmetry_file(symmetry_path)
    site_symmetries = read_site_symmetries(symmetry_path)
    representatives = sorted(atom for atom, mapped in atom_mapping.items() if atom == mapped)
    distance = float(amplitude[0])  # uniform magnitude (Phonopy convention)
    inverse_lattice = np.linalg.inv(lattice)

    entries = []
    for atom in representatives:
        rotations = np.array(site_symmetries[atom]["rotations"]) if atom in site_symmetries else np.eye(3)[None]
        for direction in _get_fractional_displacement(rotations, DIRECTIONS_DIAG):
            vector = direction @ lattice
            vector = vector / np.linalg.norm(vector) * distance
            entries.append((atom, vector))
            if is_plusminus is True or (is_plusminus == "auto" and _need_minus_fractional(direction, rotations)):
                entries.append((atom, -vector))
    return [(atom, vector @ inverse_lattice) for atom, vector in entries]


# Generate displaced coordinates and make displaced structures
def generate_displacements(
    mode: str,
    structure: Structure,
    contcar_path: str | Path,
    symmetry_path: str | Path,
    amplitude: np.ndarray,
) -> list[Displacement]:
    if mode == "minimal":
        entries = _minimal_entries(structure, symmetry_path, amplitude)
    elif mode == "fractional":
        entries = _fractional_entries(structure, symmetry_path, amplitude)
    else:
        # Determine which set of atoms to generate replacements for
        if mode == "full":
            atoms = range(1, structure.natoms + 1)
        elif mode == "atoms":
            atoms = _representatives(symmetry_path)
        else:
            raise ValueError(f"Unknown displacement mode: {mode}")
        fractional_vectors = _cartesian_vectors(amplitude) @ np.linalg.inv(structure.lattice)
        entries = [
            (atom_index, vector)
            for atom_index in atoms
            for vector in fractional_vectors
        ]

    counters: dict[int, int] = defaultdict(int)
    displacements = []
    for atom_index, vector in entries:
        suffix_index = counters[atom_index]
        if suffix_index >= len(ascii_lowercase):
            raise ValueError(f"Too many displacements for atom {atom_index}")
        displacements.append(Displacement(atom_index, vector, ascii_lowercase[suffix_index]))
        counters[atom_index] += 1
    return displacements


def _header(mode: str) -> str:
    if mode == "minimal":
        return "displacements for VASP. System (Phonopy site-symmetry-minimal set)"
    if mode == "fractional":
        return "displacements for VASP. System (Phonopy fractional-basis minimal set)"
    if mode == "atoms":
        return "displacements for VASP. Symmetry-inequivalent atoms"
    return "displacements for VASP. System"

# Write displacements to files
def write_atomic_displacements(path: str | Path, displacements: list[Displacement], amplitude: np.ndarray) -> None:
    atoms = list(dict.fromkeys(displacement.atom_index for displacement in displacements))
    with open(path, "w") as output:
        output.write(
            "Atomic displacements are "
            + " ".join(f"{value:6.3f}" for value in amplitude)
            + " Angstrom in Cartesian coordinates along x, y, and z directions\n"
        )
        output.write(f"{len(atoms):6d} number of atoms with displacements, and atom symbols and indices shown below\n")
        for atom_index in atoms:
            output.write(f"atom{atom_index:<5d}{atom_index:6d}\n")
        output.write(f"\n{len(displacements):6d} number of displacements in fractional coordinates\n")
        for displacement in displacements:
            vector = displacement.fractional_vector
            output.write(
                f"atom{displacement.atom_index:<5d}{displacement.atom_index:6d}"
                f"{vector[0]:18.12f}{vector[1]:18.12f}{vector[2]:18.12f}\n"
            )


def write_calculation_directories(displacements: list[Displacement], structure: Structure) -> None:
    for displacement in displacements:
        filename = f"pos_{displacement.name}"
        positions = structure.fractional_positions.copy()
        positions[displacement.atom_index - 1] += displacement.fractional_vector
        write_structure(filename, structure, positions)
        directory = f"ra_{filename}"
        os.makedirs(directory, exist_ok=True)
        write_structure(Path(directory) / "POSCAR", structure, positions)

# Main command chain for creating the displacement files
def run_displacements(
    mode: str = "full",
    contcar_path: str | Path = "CONTCAR",
    symmetry_path: str | Path = "symmetry",
    input_path: str | Path = "input",
) -> list[Displacement]:
    if mode not in {"full", "atoms", "minimal", "fractional"}:
        raise ValueError("mode must be full, atoms, minimal, or fractional")
    structure = read_structure(contcar_path)
    amplitude = _amplitude(input_path)
    displacements = generate_displacements(mode, structure, contcar_path, symmetry_path, amplitude)
    write_structure("ref_poscar.vasp", structure)
    write_displacements("displacements.dat", displacements, structure, _header(mode))
    write_atomic_displacements("atomic_displacement", displacements, amplitude)
    write_calculation_directories(displacements, structure)
    print(f"Generated {len(displacements)} {mode} displacement(s).")
    print("Wrote ref_poscar.vasp, displacements.dat, atomic_displacement, pos_atom*, and ra_pos_atom*/POSCAR.")
    return displacements
