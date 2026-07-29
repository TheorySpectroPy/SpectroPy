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
    label: str
    atom_index: int
    fractional_vector: np.ndarray


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
    records: list[Displacement] = []
    for line in lines[2:2 + count]:
        values = line.split()
        if len(values) < 5:
            raise ValueError(f"Invalid displacement record: {line}")
        records.append(Displacement(values[0], int(values[1]), np.array(values[2:5], dtype=float)))
    if len(records) != count:
        raise ValueError(f"Expected {count} displacement records in {path}")
    return records


def write_displacements(
    path: str | Path,
    records: list[Displacement],
    total_atoms: int,
    header: str = "displacements for VASP. System",
) -> None:
    unique_atoms = list(dict.fromkeys((record.label, record.atom_index) for record in records))
    with open(path, "w") as output:
        output.write(f"{header}\n")
        output.write(f"{len(records):6d}   Number of displacements\n")
        for record in records:
            vector = record.fractional_vector
            output.write(
                f"{record.label:<5s}{record.atom_index:6d}"
                f"{vector[0]:12.6f}{vector[1]:12.6f}{vector[2]:12.6f}\n"
            )
        output.write(f"{len(unique_atoms):6d}{total_atoms:6d}     Number of atoms in SC\n")
        for label, atom_index in unique_atoms:
            output.write(f"{label:<5s}{atom_index:6d}\n")


def suffixes(records: list[Displacement]) -> list[str]:
    counters: dict[str, int] = defaultdict(int)
    result: list[str] = []
    for record in records:
        index = counters[record.label]
        if index >= len(ascii_lowercase):
            raise ValueError(f"Too many displacements for {record.label}")
        result.append(ascii_lowercase[index])
        counters[record.label] += 1
    return result


def _amplitude(input_path: str | Path, mode: str) -> float:
    value = read_settings(input_path).displacement_amplitude
    return value if value is not None else (0.03 if mode == "minimal" else 0.01)


def _cartesian_vectors(amplitude: float) -> np.ndarray:
    return np.array([
        [amplitude, 0.0, 0.0], [-amplitude, 0.0, 0.0],
        [0.0, amplitude, 0.0], [0.0, -amplitude, 0.0],
        [0.0, 0.0, amplitude], [0.0, 0.0, -amplitude],
    ])


def _representatives(path: str | Path, total_atoms: int) -> list[int]:
    data = yaml.safe_load(Path(path).read_text()) or {}
    mapping = data.get("atom_mapping")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{path} does not contain Phonopy atom_mapping data")
    representatives = sorted({int(value) for value in mapping.values()})
    if 0 in representatives:
        representatives = [index + 1 for index in representatives]
    if not representatives or min(representatives) < 1 or max(representatives) > total_atoms:
        raise ValueError("atom_mapping indices do not match the atoms in CONTCAR")
    return representatives


def _minimal_records(structure: Structure, labels: list[str], contcar_path: str | Path, amplitude: float) -> list[Displacement]:
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
        Displacement(
            labels[int(entry["number"])],
            int(entry["number"]) + 1,
            np.asarray(entry["displacement"]) @ inverse_lattice,
        )
        for entry in phonopy.dataset["first_atoms"]
    ]


def generate_records(
    mode: str,
    structure: Structure,
    contcar_path: str | Path,
    symmetry_path: str | Path,
    amplitude: float,
) -> list[Displacement]:
    labels = atom_labels(structure.atom_symbols)
    if mode == "minimal":
        return _minimal_records(structure, labels, contcar_path, amplitude)
    if mode == "full":
        atoms = range(1, structure.natoms + 1)
    elif mode == "atoms":
        atoms = _representatives(symmetry_path, structure.natoms)
    else:
        raise ValueError(f"Unknown displacement mode: {mode}")
    fractional_vectors = _cartesian_vectors(amplitude) @ np.linalg.inv(structure.lattice)
    return [
        Displacement(labels[atom_index - 1], atom_index, vector)
        for atom_index in atoms
        for vector in fractional_vectors
    ]


def _header(mode: str) -> str:
    if mode == "minimal":
        return "displacements for VASP. System (Phonopy site-symmetry-minimal set)"
    if mode == "atoms":
        return "displacements for VASP. Symmetry-inequivalent atoms"
    return "displacements for VASP. System"


def write_atomic_displacements(path: str | Path, records: list[Displacement], amplitude: float) -> None:
    atoms = list(dict.fromkeys(record.atom_index for record in records))
    with open(path, "w") as output:
        output.write(
            f"Atomic displacements are {amplitude:6.3f} {amplitude:6.3f} "
            f"{amplitude:6.3f} Angstrom in Cartesian coordinates along x, y, and z directions\n"
        )
        output.write(f"{len(atoms):6d} number of atoms with displacements, and atom symbols and indices shown below\n")
        for atom_index in atoms:
            output.write(f"atom{atom_index:<5d}{atom_index:6d}\n")
        output.write(f"\n{len(records):6d} number of displacements in fractional coordinates\n")
        for record in records:
            vector = record.fractional_vector
            output.write(
                f"atom{record.atom_index:<5d}{record.atom_index:6d}"
                f"{vector[0]:12.6f}{vector[1]:12.6f}{vector[2]:12.6f}\n"
            )


def write_calculation_directories(records: list[Displacement], structure: Structure) -> None:
    for record, suffix in zip(records, suffixes(records)):
        filename = f"pos_atom{record.atom_index}{suffix}"
        positions = structure.fractional_positions.copy()
        positions[record.atom_index - 1] += record.fractional_vector
        write_structure(filename, structure, positions)
        directory = f"ra_{filename}"
        os.makedirs(directory, exist_ok=True)
        write_structure(Path(directory) / "POSCAR", structure, positions)


def run_displacements(
    mode: str = "full",
    contcar_path: str | Path = "CONTCAR",
    symmetry_path: str | Path = "symmetry",
    input_path: str | Path = "input",
) -> list[Displacement]:
    if mode not in {"full", "atoms", "minimal"}:
        raise ValueError("mode must be full, atoms, or minimal")
    structure = read_structure(contcar_path)
    amplitude = _amplitude(input_path, mode)
    records = generate_records(mode, structure, contcar_path, symmetry_path, amplitude)
    write_structure("ref_poscar.vasp", structure)
    write_displacements("displacements.dat", records, structure.natoms, _header(mode))
    write_atomic_displacements("atomic_displacement", records, amplitude)
    write_calculation_directories(records, structure)
    print(f"Generated {len(records)} {mode} displacement(s).")
    print("Wrote ref_poscar.vasp, displacements.dat, atomic_displacement, pos_atom*, and ra_pos_atom*/POSCAR.")
    return records
