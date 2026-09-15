from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from string import ascii_lowercase
import os
import numpy as np
import yaml

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
                f"{vector[0]:12.6f}{vector[1]:12.6f}{vector[2]:12.6f}\n"
            )
        output.write(f"{len(unique_atoms):6d}{structure.natoms:6d}     Number of atoms in SC\n")
        for atom_index in unique_atoms:
            label = labels[atom_index - 1]
            output.write(f"{label:<5s}{atom_index:6d}\n")


def _amplitude(input_path: str | Path) -> float:
    value = read_settings(input_path).displacement_amplitude
    return value if value is not None else 0.03


def _cartesian_vectors(amplitude: float) -> np.ndarray:
    return np.array([
        [amplitude, 0.0, 0.0], [-amplitude, 0.0, 0.0],
        [0.0, amplitude, 0.0], [0.0, -amplitude, 0.0],
        [0.0, 0.0, amplitude], [0.0, 0.0, -amplitude],
    ])


def _representatives(path: str | Path) -> list[int]:
    data = yaml.safe_load(Path(path).read_text()) or {}
    mapping = data.get("atom_mapping")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{path} does not contain Phonopy atom_mapping data")
    representatives = sorted({int(value) for value in mapping.values()})
    return representatives


def _minimal_entries(structure: Structure, contcar_path: str | Path, amplitude: float) -> list[tuple[int, np.ndarray]]:
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp

    phonopy = Phonopy(
        read_vasp(str(contcar_path)),
        supercell_matrix=np.eye(3, dtype=int),
        primitive_matrix=np.eye(3),
    )
    phonopy.generate_displacements(distance=amplitude, is_diagonal=False)
    inverse_lattice = np.linalg.inv(structure.lattice)
    return [
        (
            int(entry["number"]) + 1,
            np.asarray(entry["displacement"]) @ inverse_lattice,
        )
        for entry in phonopy.dataset["first_atoms"]
    ]

# Generate displaced coordinates and make displaced structures
def generate_displacements(
    mode: str,
    structure: Structure,
    contcar_path: str | Path,
    symmetry_path: str | Path,
    amplitude: float,
) -> list[Displacement]:
    if mode == "minimal":
        entries = _minimal_entries(structure, contcar_path, amplitude)
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
    if mode == "atoms":
        return "displacements for VASP. Symmetry-inequivalent atoms"
    return "displacements for VASP. System"

# Write displacements to files
def write_atomic_displacements(path: str | Path, displacements: list[Displacement], amplitude: float) -> None:
    atoms = list(dict.fromkeys(displacement.atom_index for displacement in displacements))
    with open(path, "w") as output:
        output.write(
            f"Atomic displacements are {amplitude:6.3f} {amplitude:6.3f} "
            f"{amplitude:6.3f} Angstrom in Cartesian coordinates along x, y, and z directions\n"
        )
        output.write(f"{len(atoms):6d} number of atoms with displacements, and atom symbols and indices shown below\n")
        for atom_index in atoms:
            output.write(f"atom{atom_index:<5d}{atom_index:6d}\n")
        output.write(f"\n{len(displacements):6d} number of displacements in fractional coordinates\n")
        for displacement in displacements:
            vector = displacement.fractional_vector
            output.write(
                f"atom{displacement.atom_index:<5d}{displacement.atom_index:6d}"
                f"{vector[0]:12.6f}{vector[1]:12.6f}{vector[2]:12.6f}\n"
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
    if mode not in {"full", "atoms", "minimal"}:
        raise ValueError("mode must be full, atoms, or minimal")
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
