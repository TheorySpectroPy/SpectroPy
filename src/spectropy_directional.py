"""Measured-direction raw output, with explicit physical units."""
import json
from pathlib import Path
import numpy as np
from reconstruct_dielectric_derivatives import (
    fractional_finite_differences, reconstruct_fractional_tensor,
    expand_to_equivalent_atoms, frac_to_cart_rotation)

SCHEMA = 'spectropy.directional_derivatives.v1'


def directional_document(energy, structure, measurements, sites, mapping, transfers,
                         scheme='central', equilibrium=None, bases=None):
    lattice = np.asarray(structure.lattice)
    factor = abs(np.linalg.det(lattice)) / (4*np.pi)
    samples = []
    eq = None if equilibrium is None else equilibrium[0] + 1j*equilibrium[1]
    for atom, values in measurements.items():
        rotations = sites.get(atom, {'rotations': np.eye(3)[None]})['rotations']
        measured = [(v, real + 1j*imag) for v, real, imag in values]
        for index, (unit, response, length, generated) in enumerate(
                fractional_finite_differences(measured, rotations, lattice,
                                             difference_scheme=scheme, equilibrium=eq), 1):
            cart = unit @ lattice
            scale = np.linalg.norm(cart)
            response = response * factor / scale
            samples.append(dict(atom=int(atom), direction_index=index,
                unit_cartesian=(cart/scale).tolist(), unit_fractional=unit.tolist(),
                displacement_length_angstrom=float(length*scale),
                negative_from_symmetry=generated,
                calculated_signs=1 if generated or scheme == 'forward' else 2,
                real=response.real.tolist(), imaginary=response.imag.tolist()))
    if not samples:
        raise ValueError('No measured directional responses')
    return dict(schema=SCHEMA, energy_eV=float(energy),
        response_units='angstrom^2', dielectric_basis='cartesian',
        displacement_basis='measured unit directions',
        volume_factor_angstrom3=float(factor), difference_scheme=scheme,
        lattice_angstrom=lattice.tolist(), natoms=int(structure.natoms),
        samples=samples,
        reconstruction=dict(site_rotations_fractional={str(a):np.asarray(sites.get(a,
            {'rotations': np.eye(3)[None]})['rotations']).tolist() for a in measurements},
            atom_mapping={str(a):int(r) for a,r in mapping.items()},
            transfers_fractional=[dict(representative=int(r), atom=int(a), rotation=np.asarray(t).tolist())
                                  for (r,a),t in transfers.items()], bases=bases))


def save_directional(path, document):
    Path(path).write_text(json.dumps(document, indent=2, allow_nan=False) + '\n')


def reconstruct_directional_file(path, total_atoms):
    """Reconstruct in memory for spectrum consumers; never save a raw xyz file."""
    doc = json.loads(Path(path).read_text())
    if doc.get('schema') != SCHEMA or doc.get('natoms') != total_atoms:
        raise ValueError('Invalid directional derivative schema or atom count')
    if doc.get('response_units') != 'angstrom^2':
        raise ValueError('Unexpected directional response units')
    lattice = np.asarray(doc['lattice_angstrom'], dtype=float)
    meta = doc['reconstruction']
    by_atom = {}
    seen = set()
    for sample in doc['samples']:
        atom = sample['atom']
        key = (atom, sample['direction_index'])
        if key in seen or not 1 <= atom <= total_atoms:
            raise ValueError('Duplicate or invalid directional sample')
        seen.add(key)
        unit = np.asarray(sample['unit_fractional'], dtype=float)
        cart = unit @ lattice
        response = np.asarray(sample['real']) + 1j*np.asarray(sample['imaginary'])
        if (response.shape != (3,3) or not np.isfinite(response).all()
                or not np.allclose(cart/np.linalg.norm(cart), sample['unit_cartesian'], atol=1e-10)):
            raise ValueError('Invalid directional response or direction')
        by_atom.setdefault(atom, []).append((unit, response*np.linalg.norm(cart)))
    representatives = {}
    for atom, slices in by_atom.items():
        bases = meta.get('bases')
        basis = None if bases is None else bases['atoms'][str(atom)]
        rotations = np.asarray(meta['site_rotations_fractional'][str(atom)])
        # The established solver uses real arrays; solve both parts explicitly.
        real = reconstruct_fractional_tensor([(u,t.real) for u,t in slices], rotations, lattice, basis)
        imag = reconstruct_fractional_tensor([(u,t.imag) for u,t in slices], rotations, lattice, basis)
        representatives[atom] = real + 1j*imag
    mapping = {int(a):r for a,r in meta['atom_mapping'].items()}
    transfers = {(v['representative'],v['atom']):frac_to_cart_rotation(np.asarray(v['rotation']),lattice)
                 for v in meta['transfers_fractional']}
    tensors = expand_to_equivalent_atoms(representatives,mapping,transfers)
    if set(tensors) != set(range(1,total_atoms+1)):
        raise ValueError('Directional reconstruction does not cover all atoms')
    result = np.array([tensors[a].transpose(2,0,1) for a in range(1,total_atoms+1)])
    return result.real, result.imag
