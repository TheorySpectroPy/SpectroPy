import os

import numpy as np

from .dielectric import (collect_vasprun_files,
    read_displaced_dielectrics_many, read_dielectrics_xml)
from .symmetry import read_symmetry
from .reconstruction import (prepare_reconstruction, calculate_directional_derivatives,
                             reconstruct_atomic_field, read_fractional_bases, project_crystal_field)
from .config import DEFAULT_INPUT, resolve_settings
from .displacements import read_displacement_reduction, read_displacements
from .structure import read_structure
from .paths import displacement_root
from .directional import directional_document, save_directional


def _write_dielectrics(path, laser_energy_ev, measurements, equilibrium=None):
    """Export sampled signed-displacement responses, without projection/scaling."""
    with open(path, "w") as output:
        output.write(f"# Dielectric samples for requested photon energy {laser_energy_ev:.16e} eV\n")
        output.write("# Nearest XML energy-grid samples; dimensionless Cartesian x,y,z matrices.\n")
        output.write("# No denoising or volume factor is applied. All nine reader components are retained.\n")

        def block(label, direction, real, imaginary):
            output.write(f"\n{label}\n")
            output.write("# signed_displacement_cartesian_angstrom "
                         + " ".join(f"{value:.16e}" for value in direction) + "\n")
            for name, matrix in (("Real", real), ("Imaginary", imaginary)):
                output.write(f"# {name}\n")
                for row in matrix:
                    output.write(" ".join(f"{value:.16e}" for value in row) + "\n")

        if equilibrium is not None:
            block("Equilibrium", np.zeros(3), *equilibrium)
        for atom, samples in measurements.items():
            for index, (direction, real, imaginary) in enumerate(samples, 1):
                block(f"Atom {atom} Sample {index}", direction, real, imaginary)


def _write_derivatives(path, laser_energy_ev, lattice, positions, real_tensors, imaginary_tensors):
    factor = abs(np.linalg.det(lattice)) / (4 * np.pi)
    with open(path, "w") as output:
        output.write(f"! epsilon derivatives calculated for laser frequency {laser_energy_ev:.4f} eV\n")
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
        name = f"raw_directional_derivative_{energy:.2f}"
        save_directional(os.path.join(directory, name), directional)
        print(f"Wrote {name}")


def iter_atomic_derivatives(settings, structure):
    """Read inputs once, then yield one complex Cartesian field per laser energy.

    Numerical order is explicit: finite differences -> representative solve ->
    atom transfer (or full-field projection). No tensor files are written here.
    """
    # Generation settings do not select a derivative formula. Both importers
    # feed the same measured-vector classification and reconstruction path.
    records = read_displacements()
    zero_records = [d for d in records if d.is_equilibrium]
    if len(zero_records) > 1:
        raise ValueError('Multiple equilibrium displacement records')
    displacements = [d for d in records if not d.is_equilibrium]
    equilibrium_path = settings.equilibrium_file or displacement_root() / 'ra_pos_equilibrium/vasprun.xml'
    automatic_reference = settings.equilibrium_file is None
    equilibrium_available = (bool(zero_records) or settings.equilibrium_file is not None or
                             (displacement_root() / 'ra_pos_equilibrium').exists())
    reduction = read_displacement_reduction()
    symmetry = read_symmetry(structure) if reduction != 'full' or settings.denoise else None
    sites = symmetry.sites if reduction != 'full' else {}
    mapping = symmetry.atom_mapping if reduction != 'full' else {a: a for a in range(1, structure.natoms + 1)}
    transfers = symmetry.transfers if reduction != 'full' else {}
    bases = None
    if reduction == 'minimal' and settings.saved_basis and os.path.exists('fractional_basis.json'):
        bases = read_fractional_bases('fractional_basis.json', displacements, structure.lattice, sites)
    plan = prepare_reconstruction(displacements, structure.lattice, sites, mapping, transfers,
                                  'auto', settings.denoise and reduction != 'full', bases,
                                  equilibrium_available=equilibrium_available)
    if plan.needs_equilibrium and not os.path.isfile(equilibrium_path):
        raise ValueError(f'Unpaired displacements require a completed equilibrium vasprun.xml: {equilibrium_path}')
    if plan.needs_equilibrium and automatic_reference:
        reference = read_structure(displacement_root() / 'ra_pos_equilibrium/POSCAR')
        if (reference.atom_symbols != structure.atom_symbols or
                not np.allclose(reference.lattice, structure.lattice, atol=1e-10, rtol=0)):
            raise ValueError('Equilibrium reference atom order or lattice differs from displacement reference')
        delta = reference.fractional_positions - structure.fractional_positions
        if np.max(np.linalg.norm((delta - np.rint(delta)) @ structure.lattice, axis=1)) > 1e-6:
            raise ValueError('Equilibrium reference is not the undisplaced structure')
    for atom, pairs in plan.pairs.items():
        for index, (_, _, _, negative, rotation) in enumerate(pairs, 1):
            method = ('central (measured +/-)' if negative is not None else
                      'central (symmetry-supplied -)' if rotation is not None else
                      'forward (measured signed displacement and equilibrium)')
            print(f'Atom {atom}, direction {index}: {method}')

    collect_vasprun_files(displacements)
    energies = list(settings.laser_energies)
    equilibria = (read_dielectrics_xml(equilibrium_path, energies)
                  if plan.needs_equilibrium else [None] * len(energies))
    measured_by_energy = read_displaced_dielectrics_many(displacements, structure.lattice, energies)
    for energy, equilibrium, measurements in zip(energies, equilibria, measured_by_energy):
        complex_measurements = {atom: [(direction, real + 1j * imag) for direction, real, imag in values]
                                for atom, values in measurements.items()}
        reference = None if equilibrium is None else equilibrium[0] + 1j * equilibrium[1]
        directional = calculate_directional_derivatives(complex_measurements, plan, reference)
        field = reconstruct_atomic_field(directional, plan)
        if settings.denoise and reduction == 'full':
            field = project_crystal_field(field, symmetry.operations)
        yield energy, field, directional, measurements, equilibrium


def write_atomic_outputs(energy, structure, field, directional, measurements, equilibrium, settings):
    """Convert complex arrays only at the text-output boundary."""
    parts = ({a + 1: tensor.real for a, tensor in enumerate(field)},
             {a + 1: tensor.imag for a, tensor in enumerate(field)})
    debug = None
    if settings.raw_directional_derivatives:
        debug = directional_document(energy, structure, measurements, {},
                                     'auto', equilibrium, directional=directional)
    write_derivative_outputs(energy, structure, parts, directional=debug)
    name = f'dielectric_tensor_{energy:.2f}'
    _write_dielectrics(name, energy, measurements, equilibrium)
    print(f'Wrote {name}')


def run_derivatives(input_path=DEFAULT_INPUT, use_saved_basis=None,
                    equilibrium_vasprun=None, *, settings=None):
    """Python entry point for atomic exports; the CLI uses run_tensors."""
    settings = resolve_settings(settings, input_path, saved_basis=use_saved_basis,
                                equilibrium_file=equilibrium_vasprun)
    structure = read_structure('CONTCAR')
    for energy, field, directional, measurements, equilibrium in iter_atomic_derivatives(settings, structure):
        if settings.per_atom_output:
            write_atomic_outputs(energy, structure, field, directional, measurements, equilibrium, settings)
