import os
import sys
import csv

import numpy as np

from calculate_dielectric_derivatives import (collect_vasprun_files,
    read_displaced_dielectrics_many, read_dielectrics_xml)
from process_symmetry import compute_mapping_matrices, read_site_symmetries, read_symmetry_file
from reconstruct_dielectric_derivatives import (calculate_derivatives, reconstruct_derivatives,
    reconstruct_fractional_derivatives, read_fractional_bases)
from spectropy_config import read_settings
from spectropy_displacements import read_displacement_mode, read_displacements
from spectropy_structure import read_structure
from spectropy_paths import displacement_root
from spectropy_directional import directional_document, save_directional, directional_residuals
from reconstruct_dielectric_derivatives import project_crystal_tensor_field


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
                    output.write(" ".join(f"{value:.16e}" for value in row) + "\n")


def write_derivative_outputs(energy, structure, raw, corrected=None, representatives=None, directory=".", *, directional=None):
    """Save raw measured-direction responses and the active Cartesian field.

    Raw Cartesian reconstruction is internal only. Residuals compare measured
    responses with corrected Cartesian tensors contracted onto the same seeds.
    """
    if directional is None:
        raise ValueError('Measured directional responses are required; Cartesian tensors cannot substitute for raw observations')
    def write(name, tensors):
        _write_derivatives(os.path.join(directory, name), energy, structure.lattice,
                           structure.fractional_positions, *tensors)
        print(f"Wrote {name}")
    save_directional(os.path.join(directory, f"raw_directional_derivative_{energy:.2f}.json"), directional)
    write(f"epsilon_derivative_{energy:.2f}", raw if corrected is None else corrected)
    if corrected is None:
        return []
    write(f"corrected_epsilon_derivative_{energy:.2f}", corrected)
    rows, summary = directional_residuals(directional, corrected)
    for row in rows:
        row['species'] = structure.atom_symbols[row['atom']-1]
        print(f"Atom {row['atom']}, measured direction {row['direction_index']}: "
              f"back-contracted residual = {row['residual']:.10e} (1/angstrom)")
    path = os.path.join(directory, f"site_projection_residuals_{energy:.2f}.csv")
    with open(path, 'w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save_directional(os.path.join(directory, f"directional_residual_summary_{energy:.2f}.json"), summary)
    print(f"Measured-direction relative norm-sum residual: {summary['relative_norm_sum']}")
    return rows


def run_derivatives(input_path="input", use_saved_basis=False,
                    difference_scheme=None, equilibrium_vasprun=None):
    settings = read_settings(input_path)
    site_projection = getattr(settings, 'site_symmetry_projection', False)
    dielectric_source = settings.dielectric_source
    dielectric_filename = 'vasprun.xml'
    difference_scheme = settings.difference_scheme if difference_scheme is None else difference_scheme
    automatic_reference = difference_scheme == 'forward' and equilibrium_vasprun is None
    if difference_scheme == 'forward' and equilibrium_vasprun is None:
        equilibrium_vasprun = str(displacement_root() / 'ra_pos_equilibrium' / dielectric_filename)
    if difference_scheme not in {"central", "forward"}:
        raise ValueError("difference scheme must be central or forward")
    if difference_scheme == "forward" and not os.path.isfile(equilibrium_vasprun):
        raise ValueError(f'Forward differences require the completed ra_pos_equilibrium/{dielectric_filename} calculation; run generated directories first')
    for path in ("displacements.dat", "CONTCAR"):
        if not os.path.exists(path):
            print(f"***** {path} not found *****")
            sys.exit(1)
    structure = read_structure("CONTCAR")
    if automatic_reference:
        reference = read_structure(displacement_root() / 'ra_pos_equilibrium/POSCAR')
        if (reference.atom_symbols != structure.atom_symbols or
            not np.allclose(reference.lattice, structure.lattice, atol=1e-10, rtol=0)):
            raise ValueError('Equilibrium reference atom order or lattice differs from displacement reference')
        delta = reference.fractional_positions - structure.fractional_positions
        if np.max(np.linalg.norm((delta - np.rint(delta)) @ structure.lattice, axis=1)) > 1e-6:
            raise ValueError('Equilibrium reference is not the undisplaced structure')
    displacements = read_displacements("displacements.dat")
    mode = read_displacement_mode("displacements.dat")
    if mode in {"atoms", "minimal"} or site_projection:
        if not os.path.exists("symmetry"):
            print("***** symmetry not found *****")
            sys.exit(1)
        rotations, translations, atom_mapping = read_symmetry_file("symmetry")
        _, _, mapping_matrices = compute_mapping_matrices(structure.fractional_positions, rotations, translations, atom_mapping)
        site_symmetries = read_site_symmetries("symmetry")
    bases = None
    if mode == "minimal":
        if use_saved_basis and os.path.exists("fractional_basis.json"):
            bases = read_fractional_bases("fractional_basis.json", displacements, structure.lattice, site_symmetries)
        else:
            # Reconstruct the independent directions explicitly from the
            # measured seeds and their site-symmetry images. This path also
            # supports imported displacement sets with no generated lookup.
            print("Reconstructing independent directions from the displacement set (no lookup required).")
    print(f"Collecting {dielectric_filename} files (dielectric_source: {dielectric_source})...")
    collect_vasprun_files(displacements, source=dielectric_source)
    energies = list(settings.laser_energies)
    print("Laser energies (eV): " + ", ".join(f"{energy:.2f}" for energy in energies))
    equilibria = (read_dielectrics_xml(equilibrium_vasprun, energies)
                  if difference_scheme == "forward" else [None] * len(energies))
    measured_by_energy = read_displaced_dielectrics_many(
        displacements, structure.lattice, energies, source=dielectric_source)
    for energy, equilibrium, measurements in zip(energies, equilibria, measured_by_energy):
        _write_dielectric_log(f"dielectric_tensor_{energy:.2f}", energy, measurements, structure.fractional_positions)
        def reconstruct(projected):
            sites = site_symmetries if mode != "full" else {}
            mapping = atom_mapping if mode != "full" else {atom: atom for atom in measurements}
            transfers = mapping_matrices if mode != "full" else {}
            return reconstruct_fractional_derivatives(measurements, structure.lattice, sites, mapping, transfers,
                          difference_scheme=difference_scheme, equilibrium=equilibrium,
                          site_symmetry_projection=projected, bases=bases)
        sites = site_symmetries if mode != 'full' else {}
        mapping = atom_mapping if mode != 'full' else {a:a for a in measurements}
        transfers = mapping_matrices if mode != 'full' else {}
        directional = directional_document(energy, structure, measurements, sites,
            mapping, transfers, difference_scheme, equilibrium, bases)
        if site_projection and mode=='full':
            raw=reconstruct(False)
            field=np.array([raw[0][a]+1j*raw[1][a] for a in range(1,structure.natoms+1)])
            projected=project_crystal_tensor_field(field,structure.lattice,structure.fractional_positions,
                rotations,translations,species=structure.atom_symbols)
            corrected=({a:projected[a-1].real for a in range(1,structure.natoms+1)},
                       {a:projected[a-1].imag for a in range(1,structure.natoms+1)})
        else:
            raw = reconstruct(False) if not site_projection else None
            corrected = reconstruct(True) if site_projection else None
        write_derivative_outputs(energy, structure, raw, corrected, measurements.keys(), directional=directional)
