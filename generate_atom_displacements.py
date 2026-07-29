"""Generate +/- Cartesian displacements for symmetry-inequivalent atoms only.

This is the middle ground between the ``full`` mode (all atoms, six
displacements each) and Phonopy's ``minimal`` mode (only enough directions to
span the site-symmetry-reduced displacement space).  It mirrors the Raman
pipeline's symmetry displacement generator: retain all +/- Cartesian axes, but
perform them only for one representative of each symmetry-equivalent atom set.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import yaml


def _read_poscar(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, list[str]]:
    """Return symbols, counts, lattice vectors, and coordinate lines."""
    lines = path.read_text().splitlines(keepends=True)
    if len(lines) < 8:
        raise ValueError(f"{path} is not a complete VASP POSCAR/CONTCAR file")
    scale = float(lines[1].strip())
    lattice = np.array([list(map(float, line.split()[:3])) for line in lines[2:5]]) * scale
    try:
        counts = np.array(lines[6].split(), dtype=int)
        symbols = lines[5].split()
    except ValueError:  # VASP 4 format: no element-symbol line
        counts = np.array(lines[5].split(), dtype=int)
        symbols = [f"Elem{i + 1}" for i in range(len(counts))]
    return symbols, counts, lattice, lines


def _representatives(symmetry_path: Path, total_atoms: int) -> list[int]:
    """Read Phonopy atom representatives and return 1-based atom indices."""
    data = yaml.safe_load(symmetry_path.read_text()) or {}
    mapping = data.get("atom_mapping")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{symmetry_path} does not contain Phonopy atom_mapping data")

    representatives = sorted({int(value) for value in mapping.values()})
    # Phonopy's atom_mapping is normally 0-based. Support a 1-based mapping
    # too, since historical symmetry files have used both conventions.
    if 0 in representatives:
        representatives = [index + 1 for index in representatives]
    if not representatives or min(representatives) < 1 or max(representatives) > total_atoms:
        raise ValueError("atom_mapping indices do not match the atoms in CONTCAR")
    return representatives


def run_generate_atoms(
    contcar_path: str = "CONTCAR",
    symmetry_path: str = "symmetry",
    amplitude: float = 0.01,
) -> None:
    """Write ``displacements.dat`` for every +/- Cartesian axis of each representative."""
    contcar = Path(contcar_path)
    symmetry = Path(symmetry_path)
    if not contcar.is_file():
        raise FileNotFoundError(f"CONTCAR not found: {contcar}")
    if not symmetry.is_file():
        raise FileNotFoundError(
            f"Phonopy symmetry file not found: {symmetry}. Run `phonopy --symmetry -c CONTCAR > symmetry`."
        )
    if amplitude <= 0:
        raise ValueError("displacement amplitude must be positive")

    symbols, counts, lattice, _ = _read_poscar(contcar)
    total_atoms = int(counts.sum())
    representatives = _representatives(symmetry, total_atoms)
    inverse_lattice = np.linalg.inv(lattice)
    fractional_displacements = np.array([
        [amplitude, 0.0, 0.0], [-amplitude, 0.0, 0.0],
        [0.0, amplitude, 0.0], [0.0, -amplitude, 0.0],
        [0.0, 0.0, amplitude], [0.0, 0.0, -amplitude],
    ]) @ inverse_lattice

    atom_labels: list[str] = []
    for symbol, count in zip(symbols, counts):
        atom_labels.extend(f"{symbol}{number}" for number in range(1, int(count) + 1))

    shutil.copyfile(contcar, "ref_poscar.vasp")
    with open("displacements.dat", "w") as output:
        output.write("displacements for VASP. Symmetry-inequivalent atoms\n")
        output.write(f"{len(representatives) * 6:6d}   Number of displacements\n")
        for atom_index in representatives:
            for displacement in fractional_displacements:
                output.write(
                    f"{atom_labels[atom_index - 1]:<5s}{atom_index:6d}"
                    f"{displacement[0]:12.6f}{displacement[1]:12.6f}{displacement[2]:12.6f}\n"
                )
        output.write(f"{len(representatives):6d}{total_atoms:6d}     Number of atoms in SC\n")
        for atom_index in representatives:
            output.write(f"{atom_labels[atom_index - 1]:<5s}{atom_index:6d}\n")

    print(
        f"Generated {len(representatives) * 6} +/- Cartesian displacements for "
        f"{len(representatives)} symmetry-inequivalent atom(s) "
        f"(out of {total_atoms} total atoms)."
    )
    print("Wrote ref_poscar.vasp and displacements.dat.")


if __name__ == "__main__":  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contcar", default="CONTCAR")
    parser.add_argument("--symmetry", default="symmetry")
    parser.add_argument("--amplitude", type=float, default=0.01)
    arguments = parser.parse_args()
    run_generate_atoms(arguments.contcar, arguments.symmetry, arguments.amplitude)
