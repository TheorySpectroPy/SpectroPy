import os
import sys

import numpy as np

from calculate_dielectric_derivatives import (collect_vasprun_files,
    read_displaced_dielectrics_many, read_dielectrics_xml)
from process_symmetry import compute_mapping_matrices, read_site_symmetries, read_symmetry_file
from reconstruct_dielectric_derivatives import reconstruct_fractional_derivatives, read_fractional_bases
from spectropy_config import DEFAULT_INPUT, read_settings
from spectropy_displacements import read_displacement_mode, read_displacements
from spectropy_structure import read_structure
from spectropy_paths import displacement_root
from spectropy_directional import directional_document, save_directional
from reconstruct_dielectric_derivatives import project_crystal_tensor_field


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


def write_derivative_outputs(energy, structure, final, directory=".", *, directional=None):
    """Save one final Cartesian field and, optionally, raw debugging data."""
    name = f"epsilon_derivative_{energy:.2f}"
    _write_derivatives(os.path.join(directory, name), energy, structure.lattice,
                       structure.fractional_positions, *final)
    print(f"Wrote {name}")
    if directional is not None:
        save_directional(os.path.join(directory, f"raw_directional_derivative_{energy:.2f}.json"), directional)


def run_derivatives(input_path=DEFAULT_INPUT, use_saved_basis=None,
                    difference_scheme=None, equilibrium_vasprun=None, *, on_reconstructed=None):
    """Reconstruct one energy at a time; export atomic data only when requested."""
    settings = read_settings(input_path)
    use_saved_basis = settings.saved_basis if use_saved_basis is None else use_saved_basis
    equilibrium_vasprun = settings.equilibrium_file if equilibrium_vasprun is None else equilibrium_vasprun
    site_projection = getattr(settings, 'site_symmetry_projection', False)
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
    print(f"Collecting {dielectric_filename} files...")
    collect_vasprun_files(displacements, source=dielectric_filename)
    energies = list(settings.laser_energies)
    print("Laser energies (eV): " + ", ".join(f"{energy:.2f}" for energy in energies))
    equilibria = (read_dielectrics_xml(equilibrium_vasprun, energies)
                  if difference_scheme == "forward" else [None] * len(energies))
    measured_by_energy = read_displaced_dielectrics_many(
        displacements, structure.lattice, energies, source=dielectric_filename)
    for energy, equilibrium, measurements in zip(energies, equilibria, measured_by_energy):
        def reconstruct(projected):
            sites = site_symmetries if mode != "full" else {}
            mapping = atom_mapping if mode != "full" else {atom: atom for atom in measurements}
            transfers = mapping_matrices if mode != "full" else {}
            return reconstruct_fractional_derivatives(measurements, structure.lattice, sites, mapping, transfers,
                          difference_scheme=difference_scheme, equilibrium=equilibrium,
                          site_symmetry_projection=projected, bases=bases)
        directional = None
        if settings.raw_directional_derivatives:
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
            final=({a:projected[a-1].real for a in range(1,structure.natoms+1)},
                   {a:projected[a-1].imag for a in range(1,structure.natoms+1)})
        else:
            final = reconstruct(site_projection)
        if on_reconstructed is not None:
            on_reconstructed(energy, structure, final)
        if settings.per_atom_output:
            write_derivative_outputs(energy, structure, final, directional=directional)
