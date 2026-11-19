"""Import VASP/Angstrom single-atom Phonopy datasets without reselecting directions."""
from collections import defaultdict
from pathlib import Path
from string import ascii_lowercase

import numpy as np

from .displacements import Displacement
from .reconstruction import prepare_reconstruction, build_fractional_bases


def read_phonopy_displacements(path, structure):
    """Read type-1 records in YAML order; retain their physical vectors exactly.

    The reference must be the YAML supercell, not its primitive/unit cell.
    Phonopy's reader handles YAML's 1-based to Python's 0-based atom indices.
    Only the VASP interface is supported here; other calculator units must not
    silently be interpreted as Angstrom.
    """
    try:
        from phonopy.interface.phonopy_yaml import PhonopyYaml
    except ImportError as error:
        raise ImportError('Phonopy import requires the phonopy package: pip install "spectropy[phonopy]"') from error
    data = PhonopyYaml().read(str(path))
    if data.calculator not in (None, 'vasp'):
        raise ValueError('Phonopy import currently supports VASP/Angstrom datasets only')
    units = data.physical_units
    length_unit = (units.get('length') if isinstance(units, dict)
                   else getattr(units, 'length_unit', None))
    if length_unit is not None and str(length_unit).lower() != 'angstrom':
        raise ValueError('Phonopy displacement lengths must be in Angstrom')
    cell = data.supercell
    if cell is None:
        raise ValueError('Phonopy YAML must include the undisplaced supercell')
    if (list(cell.symbols) != structure.atom_symbols or
            not np.allclose(cell.cell, structure.lattice, atol=1e-8, rtol=0)):
        raise ValueError('CONTCAR atom order/species or lattice differs from the Phonopy supercell')
    delta = cell.scaled_positions - structure.fractional_positions
    if not np.isfinite(delta).all() or np.max(np.linalg.norm((delta - np.rint(delta)) @ structure.lattice, axis=1)) > 1e-7:
        raise ValueError('CONTCAR positions differ from the undisplaced Phonopy supercell')
    dataset = data.dataset
    if not dataset or 'first_atoms' not in dataset:
        raise ValueError('Only single-atom Phonopy displacement datasets are supported; random/multi-atom snapshots are not')
    if dataset.get('natom') != structure.natoms or not dataset['first_atoms']:
        raise ValueError('Empty or inconsistent Phonopy displacement dataset')
    inverse = np.linalg.inv(structure.lattice)
    counters, previous = defaultdict(int), defaultdict(list)
    records = []
    for entry in dataset['first_atoms']:
        index = entry['number']
        if isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)) or not 0 <= index < structure.natoms:
            raise ValueError(f'Invalid Phonopy atom index: {index}')
        atom = int(index) + 1
        vector = np.asarray(entry['displacement'], dtype=float)
        if vector.shape != (3,) or not np.isfinite(vector).all():
            raise ValueError(f'Atom {atom}: displacement must be a finite three-vector')
        length = np.linalg.norm(vector)
        if not np.isfinite(length):
            raise ValueError(f'Atom {atom}: displacement length must be finite')
        if length == 0:
            if any(d.is_equilibrium for d in records):
                raise ValueError('Repeated equilibrium displacement')
            records.append(Displacement(atom, np.zeros(3), ''))
            continue
        for old in previous[atom]:
            old_length = np.linalg.norm(old)
            dot = np.dot(vector / length, old / old_length)
            if abs(dot) > 1 - 1e-10:
                if dot > 0:
                    raise ValueError(f'Atom {atom}: repeated displacement direction; import one amplitude at a time')
                if not np.isclose(length, old_length, atol=1e-10, rtol=1e-8):
                    raise ValueError(f'Atom {atom}: opposite displacements must have equal lengths')
        suffix = counters[atom]
        if suffix >= len(ascii_lowercase):
            raise ValueError(f'Too many displacements for atom {atom}')
        records.append(Displacement(atom, vector @ inverse, ascii_lowercase[suffix]))
        previous[atom].append(vector)
        counters[atom] += 1
    return records


def validate_phonopy_reconstruction(records, structure, symmetry, reduction, settings):
    """Reject incomplete directions/signs or incompatible representatives before writes."""
    if reduction == 'full':
        mapping = {a: a for a in range(1, structure.natoms + 1)}
        sites, transfers = {}, {}
    else:
        mapping, sites, transfers = symmetry.atom_mapping, symmetry.sites, symmetry.transfers
    equilibrium = any(d.is_equilibrium for d in records) or settings.equilibrium_file is not None
    records = [d for d in records if not d.is_equilibrium]
    missing = set(mapping.values()) - {d.atom_index for d in records}
    if missing:
        raise ValueError(f'Phonopy set is missing representative atoms {sorted(missing)}; '
                         'use matching supercell symmetry and symmetry.reduction: on for a reduced set')
    # Match the vectors that displacements.dat will actually serialize.
    rounded = [Displacement(d.atom_index, np.round(d.fractional_vector, 12), d.suffix) for d in records]
    def prepare(has_equilibrium):
        return prepare_reconstruction(rounded, structure.lattice, sites, mapping, transfers,
                                      'auto', settings.denoise and reduction != 'full',
                                      equilibrium_available=has_equilibrium)
    plan = prepare(equilibrium)
    if plan.needs_equilibrium and not equilibrium:
        # A one-sided set needs the original cell as an optical reference.
        # Once present, it is used for all unpaired measured directions.
        plan = prepare(True)
    if settings.saved_basis and reduction == 'minimal':
        build_fractional_bases(rounded, structure.lattice, sites)
    return plan


def check_import_destination(root):
    """Imports only stage new datasets; never relabel/overwrite existing calculations."""
    for name in ('displacements.dat', 'atomic_displacement',
                 'fractional_basis.json', 'phonopy_displacement_map.tsv'):
        path = Path(name)
        if path.exists() or path.is_symlink():
            raise FileExistsError(f'Refusing to overwrite {path}; import in a fresh working directory')
    # Workflow launchers and input links may already occupy this container.
    # Protect the generated names, not unrelated files we will never write.
    occupied = (root.is_symlink() or (root.exists() and not root.is_dir())
                or any(root.glob('ra_pos_*')) or any(root.glob('pos_atom*'))
                or (root / 'vasprun').exists() or (root / 'vasprun').is_symlink())
    if occupied:
        raise FileExistsError(f'Refusing to overwrite {root}; import in a fresh working directory')


def write_phonopy_mapping(path, records, lattice):
    """Map YAML record order to generated calculation names (not a solver input)."""
    with open(path, 'w') as stream:
        stream.write('phonopy_record\tatom\tdx_A\tdy_A\tdz_A\tspectropy_directory\n')
        for index, record in enumerate(records, 1):
            vector = record.fractional_vector @ lattice
            stream.write(f'{index}\t{record.atom_index}\t'
                         + '\t'.join(f'{x:.16e}' for x in vector)
                         + f'\tra_pos_{record.name}\n')
