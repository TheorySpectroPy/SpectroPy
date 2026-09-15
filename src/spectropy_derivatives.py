import os
import sys

import numpy as np

from calculate_dielectric_derivatives import collect_vasprun_files, read_displaced_dielectrics
from process_symmetry import compute_mapping_matrices, read_site_symmetries, read_symmetry_file
from reconstruct_dielectric_derivatives import calculate_derivatives, reconstruct_derivatives
from spectropy_config import read_settings
from spectropy_displacements import read_displacement_mode, read_displacements
from spectropy_structure import read_structure


def _write_dielectric_log(path, target_frequency, measurements, positions):
    with open(path, "w") as output:
        output.write(f" Atomic displacement dielectric tensors at {target_frequency:.2f} eV\n")
        for atom, values in measurements.items():
            for direction, real, imaginary in values:
                position = positions[atom - 1]
                output.write(f" atom{atom} " + " ".join(f"{value:9.6f}" for value in position) + " " + " ".join(f"{value:9.6f}" for value in direction) + "\n")
                for row in range(3):
                    output.write("".join(f" {value:12.7f} {imaginary[row, column]:12.7f}" for column, value in enumerate(real[row])) + "\n")


def _write_derivatives(path, target_frequency, lattice, positions, real_tensors, imaginary_tensors):
    factor = np.linalg.det(lattice) / (4 * np.pi)
    with open(path, "w") as output:
        output.write(f"! epsilon derivatives calculated for laser frequency {target_frequency:.4f} eV\n")
        output.write("! Unit cell matrix:\n")
        for vector in lattice.T:
            output.write(f"!   {vector[0]:21.16f} {vector[1]:21.16f} {vector[2]:21.16f}\n")
        output.write("!======================================================\n!-------------------------------------------------------\n")
        for title, tensors in ((" ! Real Part of epsilon derivative tensor per atom along x, y, and z directions", real_tensors), ("\n ! Imaginary Part of epsilon derivative tensor per atom along x, y, and z directions", imaginary_tensors)):
            output.write(title + "\n")
            for atom in sorted(tensors):
                position = positions[atom - 1]
                output.write(f"      Atom {atom} {position[0]:10.6f} {position[1]:10.6f} {position[2]:10.6f}\n")
                matrix = np.hstack(tuple(tensors[atom][:, :, axis] * factor for axis in range(3)))
                for row in matrix:
                    output.write("".join(f"{value:16.4f}" for value in row) + "\n")


def run_derivatives(input_path="input"):
    for path in ("displacements.dat", "CONTCAR"):
        if not os.path.exists(path):
            print(f"***** {path} not found *****")
            sys.exit(1)
    structure = read_structure("CONTCAR")
    displacements = read_displacements("displacements.dat")
    mode = read_displacement_mode("displacements.dat")
    if mode in {"atoms", "minimal"}:
        if not os.path.exists("symmetry"):
            print("***** symmetry not found *****")
            sys.exit(1)
        rotations, translations, atom_mapping = read_symmetry_file("symmetry")
        _, _, mapping_matrices = compute_mapping_matrices(structure.fractional_positions, rotations, translations, atom_mapping)
        site_symmetries = read_site_symmetries("symmetry")
    print("Collecting vasprun.xml files...")
    collect_vasprun_files(displacements)
    energies = read_settings(input_path).laser_energies
    print("Laser energies (eV): " + ", ".join(f"{energy:.2f}" for energy in energies))
    for energy in energies:
        measurements = read_displaced_dielectrics(displacements, structure.lattice, energy)
        _write_dielectric_log(f"dielectric_tensor_{energy:.2f}", energy, measurements, structure.fractional_positions)
        if mode == "full":
            real, imaginary = calculate_derivatives(measurements, structure.lattice)
        else:
            real, imaginary = reconstruct_derivatives(measurements, structure.lattice, site_symmetries, atom_mapping, mapping_matrices)
        _write_derivatives(f"epsilon_derivative_{energy:.2f}", energy, structure.lattice, structure.fractional_positions, real, imaginary)
        print(f"Wrote epsilon_derivative_{energy:.2f}")
