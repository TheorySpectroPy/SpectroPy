"""Read-only regression on archived fractional dielectric logs.

Run with PYTHONPATH=src python tests/check_fractional_basis_data.py RAMAN_DIR [...].
The exact displacement-file vectors are used, not rounded log directions.
No calculation input or output is modified. JSON round-trip uses a temp dir.
"""
from collections import defaultdict
import hashlib
from pathlib import Path
import sys
import tempfile

import numpy as np

from process_symmetry import read_site_symmetries, read_symmetry_file, compute_mapping_matrices
from reconstruct_dielectric_derivatives import (build_fractional_bases, write_fractional_bases,
    read_fractional_bases, reconstruct_fractional_derivatives)
from spectropy_displacements import read_displacements
from spectropy_structure import read_structure
from spectropy_phonons import read_gamma_modes

_mode_cache = {}


def check(directory):
    structure = read_structure(directory / "CONTCAR")
    displacements = read_displacements(directory / "displacements.dat")
    site = read_site_symmetries(directory / "symmetry")
    rotations, translations, mapping = read_symmetry_file(directory / "symmetry")
    _, _, matrices = compute_mapping_matrices(structure.fractional_positions, rotations, translations, mapping)
    bases = build_fractional_bases(displacements, structure.lattice, site)
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "fractional_basis.json"
        write_fractional_bases(path, bases)
        bases = read_fractional_bases(path, displacements, structure.lattice, site)
    band_path = directory / "band.yaml"
    checksum = hashlib.sha256(band_path.read_bytes()).hexdigest()
    if checksum not in _mode_cache:
        _mode_cache[checksum] = read_gamma_modes(band_path)
    modes = _mode_cache[checksum]
    for log in sorted(directory.glob("dielectric_tensor_*")):
        lines = log.read_text().splitlines()[1:]
        assert len(lines) == 4 * len(displacements), log
        measurements = defaultdict(list)
        for index, displacement in enumerate(displacements):
            assert int(lines[4 * index].split()[0][4:]) == displacement.atom_index
            rows = np.array([list(map(float, line.split())) for line in lines[4 * index + 1:4 * index + 4]])
            measurements[displacement.atom_index].append((displacement.fractional_vector @ structure.lattice, rows[:, ::2], rows[:, 1::2]))
        args = (measurements, structure.lattice, site, mapping, matrices)
        legacy = reconstruct_fractional_derivatives(*args)
        saved = reconstruct_fractional_derivatives(*args, bases=bases)
        for old, new in zip(legacy, saved):
            for atom in old:
                np.testing.assert_array_equal(old[atom], new[atom])
        old_field = np.array([legacy[0][a] + 1j * legacy[1][a] for a in sorted(mapping)])
        new_field = np.array([saved[0][a] + 1j * saved[1][a] for a in sorted(mapping)])
        old_modes = np.einsum("jabk,mjk->mab", old_field, modes.eigendisplacements)
        new_modes = np.einsum("jabk,mjk->mab", new_field, modes.eigendisplacements)
        np.testing.assert_array_equal(old_modes, new_modes)
        np.testing.assert_array_equal(np.abs(old_modes[:, 0, 0]) ** 2, np.abs(new_modes[:, 0, 0]) ** 2)
        print(f"{directory.parent.name} {log.name}: {len(mapping)} atomic tensors, {len(old_modes)} complex mode tensors identical; max difference=0")


if __name__ == "__main__":
    for argument in sys.argv[1:]:
        check(Path(argument))
