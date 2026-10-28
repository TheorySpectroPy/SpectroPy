"""Exact, orbit-wise constraints on Cartesian atomic dielectric derivatives.

This diagnostic derives forms only: it never changes numerical derivatives.
The supplied operations must be the complete unitary symmetry group of the
electronic state, represented by (atom permutation, Cartesian rotation).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def _zero(matrix):
    # Covariance differences are usually polynomial identities in the free
    # parameters. Exact expansion avoids expensive general assumption queries;
    # retain simplify as a fallback for expressions needing other identities.
    return all(sp.expand(x) == 0 or sp.simplify(x) == 0 for x in matrix)


def exact_rotation(rotation, tolerance=1e-9):
    """Recover small crystallographic radicals, rejecting non-exact isometries.

    Supply symbolic strings for orientations not recoverable with these radicals.
    No silently approximate matrix is used in the symbolic linear equations.
    """
    raw = sp.Matrix(rotation)
    if raw.shape != (3, 3):
        raise ValueError('Cartesian rotations must be 3 by 3')
    exact = raw.applyfunc(lambda x: sp.nsimplify(
        x, [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)], tolerance=tolerance, full=False))
    if not _zero(exact.T * exact - sp.eye(3)):
        raise ValueError('Rotation cannot be recovered as an exact orthogonal matrix; '
                         'supply exact symbolic Cartesian entries')
    if max(abs(complex(sp.N(x))) for x in exact - raw) > tolerance:
        raise ValueError('Exact rotation differs from input beyond tolerance')
    return exact


def tensor_action(rotation):
    """27-vector order: (displacement i, dielectric alpha, dielectric beta)."""
    return sp.kronecker_product(rotation, rotation, rotation)


def derive_constraints(natoms, operations, symmetric_dielectric=False):
    """Return exact allowed forms, including shared parameters across atom orbits.

    Parameters are unconstrained scalars (complex at finite photon frequency).
    Dielectric reciprocity is optional; do not impose it blindly in magnets.
    Translation information enters through the supplied atom permutations.
    """
    if natoms < 1 or not operations:
        raise ValueError('A positive atom count and a full symmetry group are required')
    ops = []
    for permutation, rotation in operations:
        p = tuple(permutation)
        if sorted(p) != list(range(natoms)):
            raise ValueError('Each atom permutation must contain 0..natoms-1 exactly once')
        s = exact_rotation(rotation)
        if not any(p == q and _zero(s - t) for q, t in ops):
            ops.append((p, s))
    # Identity and closure are essential for orbit/site decomposition. Distinct
    # Seitz operations with identical Gamma-point actions may be deduplicated.
    if not any(p == tuple(range(natoms)) and _zero(s - sp.eye(3)) for p, s in ops):
        raise ValueError('The supplied group is missing its identity')
    keys = {(p, tuple(sp.simplify(x) for x in s)) for p, s in ops}
    for p, s in ops:
        for q, t in ops:
            product = (tuple(p[q[i]] for i in range(natoms)),
                       tuple(sp.simplify(x) for x in s * t))
            if product not in keys:
                raise ValueError('Operations are not a complete closed group')
    actions = [tensor_action(s) for _, s in ops]
    forms, orbits = {}, []
    for representative in range(natoms):
        if representative in forms:
            continue
        orbit = sorted({p[representative] for p, _ in ops})
        site = [i for i, (p, _) in enumerate(ops) if p[representative] == representative]
        equations = [actions[i] - sp.eye(27) for i in site]
        if symmetric_dielectric:
            reciprocity = sp.zeros(9, 27)
            row = 0
            for i in range(3):
                for a, b in ((0, 1), (0, 2), (1, 2)):
                    reciprocity[row, 9*i + 3*a + b] = 1
                    reciprocity[row, 9*i + 3*b + a] = -1
                    row += 1
            equations.append(reciprocity)
        reduced, pivots = sp.Matrix.vstack(*equations).rref()
        free = [i for i in range(27) if i not in pivots]
        parameters = [sp.Symbol(f'p{representative+1}_{j+1}') for j in range(len(free))]
        vector = sp.zeros(27, 1)
        for column, symbol in zip(free, parameters):
            vector[column] = symbol
        for row, column in enumerate(pivots):
            vector[column] = sp.simplify(-sum(reduced[row, j]*vector[j] for j in free))
        transfers = {}
        for atom in orbit:
            index = next(i for i, (p, _) in enumerate(ops) if p[representative] == atom)
            forms[atom] = (actions[index] * vector).applyfunc(sp.simplify)
            transfers[atom] = index
        # Independent all-operation check, including every alternative atom map.
        for atom in orbit:
            for (p, _), action in zip(ops, actions):
                if not _zero(forms[p[atom]] - action * forms[atom]):
                    raise AssertionError('Derived forms fail full covariance verification')
        orbits.append(dict(representative=representative, atoms=orbit,
                           site_operations=site, parameters=parameters,
                           free_components=free, transfers=transfers))
    return dict(forms=forms, orbits=orbits, operations=ops,
                symmetric_dielectric=symmetric_dielectric,
                parameter_count=sum(len(o['parameters']) for o in orbits))


def operations_from_poscar(path, symprec=1e-5):
    """Structural unitary operations; caller must ensure electronic-state symmetry."""
    import spglib
    from spectropy_structure import read_structure
    structure = read_structure(path)
    species = np.repeat(np.arange(len(structure.counts)), structure.counts) + 1
    data = spglib.get_symmetry((structure.lattice, structure.fractional_positions, species),
                              symprec=symprec)
    if data is None:
        raise ValueError('spglib did not find structural symmetry')
    operations = []
    for w, t in zip(data['rotations'], data['translations']):
        moved = structure.fractional_positions @ w.T + t
        delta = structure.fractional_positions[None] - moved[:, None]
        delta -= np.rint(delta)
        distances = np.linalg.norm(delta @ structure.lattice, axis=2)
        distances[species[:, None] != species[None, :]] = np.inf
        permutation = distances.argmin(axis=1)
        if len(set(permutation)) != structure.natoms or np.max(
                distances[np.arange(structure.natoms), permutation]) > symprec:
            raise ValueError('Ambiguous or inaccurate symmetry atom mapping')
        rotation = structure.lattice.T @ w @ np.linalg.inv(structure.lattice.T)
        operations.append((permutation.tolist(), rotation.tolist()))
    return structure, operations


def markdown_report(result, labels=None):
    """Human-readable exact matrices and complete independent equations."""
    labels = labels or [f'atom {i+1}' for i in range(len(result['forms']))]
    lines = ['# Symbolic atomic dielectric-derivative constraints', '',
             'D_i(mu)[alpha,beta] = partial epsilon_alpha,beta / partial r_i(mu).',
             'All three indices are Cartesian. Matrix rows/columns are dielectric x,y,z;',
             'the slice index i is the displacement direction. No numerical data are modified.', '',
             'Full covariance: D(g mu) = (S_g tensor S_g tensor S_g) D(mu).',
             f"Dielectric-index symmetry imposed: {result['symmetric_dielectric']}.",
             f"Independent scalar parameters across all atoms: {result['parameter_count']}.",
             'Parameters may be complex. Operations must preserve the electronic state.', '']
    for orbit in result['orbits']:
        rep = orbit['representative']
        lines += [f'## Orbit represented by {labels[rep]} (atom {rep+1})', '',
                  f"Atoms (1-based): {[a+1 for a in orbit['atoms']]}; "
                  f"site operations: {len(orbit['site_operations'])}; "
                  f"independent parameters: {len(orbit['parameters'])}.", '',
                  'The same parameters are shared by all atoms in this orbit.', '']
        for atom in orbit['atoms']:
            vector = result['forms'][atom]
            lines += [f'### {labels[atom]} (atom {atom+1})', '']
            for i, axis in enumerate('xyz'):
                matrix = sp.Matrix(3, 3, list(vector[9*i:9*i+9]))
                lines += [f'D_{axis}:', '', '```text', str(matrix), '```', '',
                          '$$' + sp.latex(matrix) + '$$', '']
        # Parameter definitions + all other components = complete independent
        # constraints, not merely a list of entries expected to vanish.
        v = result['forms'][rep]
        free = set(orbit['free_components'])
        lines += ['### Representative component equations', '', '```text']
        for j in range(27):
            i, a, b = j//9, (j%9)//3, j%3
            name = f'D_{"xyz"[i]}[{"xyz"[a]},{"xyz"[b]}]'
            lines.append(f'{name} = {v[j]}' + (' (parameter definition)' if j in free else ''))
        lines += ['```', '']
    lines += ['## Exact Cartesian operations', '',
              'Permutations below map source atom to target atom (1-based).', '']
    for index, (p, s) in enumerate(result['operations']):
        lines += [f'### Operation {index+1}: {[i+1 for i in p]}', '',
                  '```text', str(s), '```', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--poscar', type=Path)
    source.add_argument('--operations', type=Path,
                        help='JSON: natoms, operations [{permutation: zero-based list, rotation: 3x3}]')
    parser.add_argument('--symmetric-dielectric', action='store_true')
    parser.add_argument('--symprec', type=float, default=1e-5)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; choose a new path to preserve previous reports')
    if args.poscar:
        structure, operations = operations_from_poscar(args.poscar, args.symprec)
        natoms, labels = structure.natoms, structure.atom_symbols
    else:
        data = json.loads(args.operations.read_text())
        natoms, labels = data['natoms'], data.get('labels')
        operations = [(op['permutation'], op['rotation']) for op in data['operations']]
    result = derive_constraints(natoms, operations, args.symmetric_dielectric)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(markdown_report(result, labels))
    print(f"{result['parameter_count']} independent parameters; report: {args.output}")


if __name__ == '__main__':
    main()
