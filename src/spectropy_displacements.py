from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from string import ascii_lowercase
import os
import shutil
import numpy as np
import yaml

from process_symmetry import read_site_symmetries, read_symmetry_file
from reconstruct_dielectric_derivatives import build_fractional_bases, write_fractional_bases
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


def read_displacement_metadata(path: str | Path = "displacements.dat") -> tuple[str, str]:
    fields = dict(item.split("=", 1) for item in Path(path).read_text().splitlines()[0].split()[1:])
    mode, basis = fields.get("mode"), fields.get("basis")
    if mode not in {"full", "atoms", "minimal"} or basis not in {"cartesian", "fractional"}:
        raise ValueError(f"Invalid displacement mode/basis metadata in {path}")
    return mode, basis


def read_displacement_mode(path: str | Path = "displacements.dat") -> str:
    return read_displacement_metadata(path)[0]


def write_displacements(
    path: str | Path,
    displacements: list[Displacement],
    structure: Structure,
    header: str = "SpectroPy mode=full basis=cartesian",
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
    """Physical lengths along candidate basis axes; one value applies to all."""
    value = read_settings(input_path).displacement_amplitude
    return np.array(value if value is not None else (0.03, 0.03, 0.03))


def _representatives(path: str | Path) -> list[int]:
    data = yaml.safe_load(Path(path).read_text()) or {}
    mapping = data.get("atom_mapping")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{path} does not contain Phonopy atom_mapping data")
    representatives = sorted({int(value) for value in mapping.values()})
    return representatives


# ---------------------------------------------------------------------------
# Shared basis-aware minimal selection
# ---------------------------------------------------------------------------

def _frac_to_cart_rotation(rotation: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """Cartesian-basis form of a fractional-basis rotation.

    Matches reconstruct_dielectric_derivatives.frac_to_cart_rotation so the
    generator and the derivative reconstruction share one convention.
    """
    return lattice.T @ rotation @ np.linalg.inv(lattice.T)


def _need_minus(direction, rotations, tol=1e-8):
    return not any(np.allclose(rotation @ direction, -direction, atol=tol, rtol=0)
                   for rotation in rotations)


# Candidate coordinates in the selected basis, in Phonopy's directions_diag
# order. Both bases share this list and the same tolerant orbit search.
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
    """All site-symmetry images in the chosen basis."""
    return [direction @ rotation.T for rotation in rotations]


def _get_displacement_one(rotations, directions):
    """One candidate whose site-symmetry orbit spans 3D."""
    for direction in directions:
        images = _rotated(direction, rotations)
        for i in range(len(rotations)):
            for j in range(i + 1, len(rotations)):
                if abs(_det(direction, images[i], images[j])) > 1e-8:
                    return [direction]
    return None


def _get_displacement_two(rotations, directions):
    """Two candidates whose site-symmetry orbits span 3D."""
    for direction in directions:
        images = _rotated(direction, rotations)
        for i in range(len(rotations)):
            for second in directions:
                if abs(_det(direction, images[i], second)) > 1e-8:
                    return [direction, second]
    return None


def _select_directions(rotations, directions):
    """Shared one/two/three search; not a global signed-cost optimum."""
    result = _get_displacement_one(rotations, directions)
    if result is not None:
        return result
    result = _get_displacement_two(rotations, directions)
    if result is not None:
        return result
    return [directions[0], directions[1], directions[2]]


# Generate displaced coordinates and make displaced structures
def generate_displacements(
    mode: str,
    structure: Structure,
    contcar_path: str | Path,
    symmetry_path: str | Path,
    amplitude: np.ndarray,
    basis: str = "cartesian",
) -> list[Displacement]:
    if mode not in {"full", "atoms", "minimal"}:
        raise ValueError("mode must be full, atoms, or minimal")
    if basis not in {"cartesian", "fractional"}:
        raise ValueError("basis must be cartesian or fractional")
    amplitude = np.asarray(amplitude, dtype=float)
    if not np.isfinite(amplitude).all() or np.any(amplitude <= 0):
        raise ValueError("Displacement amplitudes must be positive and finite")
    if mode == "minimal" and not np.allclose(amplitude, amplitude[0]):
        raise ValueError("minimal requires a common physical displacement length")
    atoms = range(1, structure.natoms + 1) if mode == "full" else _representatives(symmetry_path)
    sites = read_site_symmetries(symmetry_path) if mode == "minimal" else {}
    inverse_lattice = np.linalg.inv(structure.lattice)
    entries = []
    for atom in atoms:
        rotations = np.asarray(sites.get(atom, {}).get("rotations", [np.eye(3)]))
        if basis == "cartesian":
            rotations = np.asarray([_frac_to_cart_rotation(r, structure.lattice) for r in rotations])
        directions = _select_directions(rotations, DIRECTIONS_DIAG) if mode == "minimal" else np.eye(3)
        for index, direction in enumerate(directions):
            vector = direction @ structure.lattice if basis == "fractional" else direction
            length = amplitude[0] if mode == "minimal" else amplitude[index]
            vector = vector / np.linalg.norm(vector) * length @ inverse_lattice
            entries.append((atom, vector))
            if mode != "minimal" or _need_minus(direction, rotations):
                entries.append((atom, -vector))

    counters: dict[int, int] = defaultdict(int)
    displacements = []
    for atom_index, vector in entries:
        suffix_index = counters[atom_index]
        if suffix_index >= len(ascii_lowercase):
            raise ValueError(f"Too many displacements for atom {atom_index}")
        displacements.append(Displacement(atom_index, vector, ascii_lowercase[suffix_index]))
        counters[atom_index] += 1
    return displacements


def _header(mode: str, basis: str) -> str:
    return f"SpectroPy mode={mode} basis={basis}"

# Write displacements to files
def write_atomic_displacements(path: str | Path, displacements: list[Displacement], amplitude: np.ndarray) -> None:
    atoms = list(dict.fromkeys(displacement.atom_index for displacement in displacements))
    with open(path, "w") as output:
        output.write(
            "Atomic displacements are "
            + " ".join(f"{value:6.3f}" for value in amplitude)
            + " Angstrom physical lengths; displacement vectors below are fractional\n"
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
    root = Path('displacements') if Path.cwd().name != 'displacements' else Path('.')
    root.mkdir(exist_ok=True)
    for displacement in displacements:
        filename = f"pos_{displacement.name}"
        positions = structure.fractional_positions.copy()
        positions[displacement.atom_index - 1] += displacement.fractional_vector
        write_structure(root / filename, structure, positions)
        directory = root / f"ra_{filename}"
        os.makedirs(directory, exist_ok=True)
        write_structure(Path(directory) / "POSCAR", structure, positions)

# Main command chain for creating the displacement files
def run_displacements(
    mode: str = "full",
    contcar_path: str | Path = "CONTCAR",
    symmetry_path: str | Path = "symmetry",
    input_path: str | Path = "input",
    basis: str = "cartesian",
) -> list[Displacement]:
    if mode not in {"full", "atoms", "minimal"}:
        raise ValueError("mode must be full, atoms, or minimal")
    structure = read_structure(contcar_path)
    amplitude = _amplitude(input_path)
    displacements = generate_displacements(mode, structure, contcar_path, symmetry_path, amplitude, basis=basis)
    write_structure("ref_poscar.vasp", structure)
    write_displacements("displacements.dat", displacements, structure, _header(mode, basis))
    if mode == "minimal":
        bases = build_fractional_bases(read_displacements("displacements.dat"), structure.lattice, read_site_symmetries(symmetry_path))
        write_fractional_bases("fractional_basis.json", bases)
        print("Saved fractional reconstruction bases to fractional_basis.json.")
    write_atomic_displacements("atomic_displacement", displacements, amplitude)
    write_calculation_directories(displacements, structure)
    # Include the undisplaced optical reference in the workflow's ra_pos_* loop,
    # but not in the nonzero displacement records used for differentiation.
    directory = (Path('displacements') if Path.cwd().name != 'displacements' else Path('.')) / 'ra_pos_equilibrium'
    directory.mkdir(exist_ok=True)
    path = directory / 'POSCAR'
    if path.exists():
        existing = read_structure(path)
        delta = existing.fractional_positions - structure.fractional_positions if existing.natoms == structure.natoms else None
        if (existing.atom_symbols != structure.atom_symbols or
            not np.allclose(existing.lattice, structure.lattice, atol=1e-10, rtol=0) or
            delta is None or not np.allclose(delta - np.rint(delta), 0, atol=1e-10, rtol=0)):
            raise ValueError('Existing equilibrium POSCAR differs from displacement reference; refusing to overwrite')
    else:
        shutil.copyfile(contcar_path, path)
    print(f'Generated {directory}/POSCAR (undisplaced optical reference).')
    print(f"Generated {len(displacements)} {mode} displacement(s).")
    print("Wrote reference metadata and displacements/{pos_atom*,ra_pos_atom*/POSCAR}.")
    return displacements
