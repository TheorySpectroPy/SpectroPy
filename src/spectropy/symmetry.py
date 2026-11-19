import numpy as np
import os
import sys
import yaml
from dataclasses import dataclass
from .structure import read_structure


def frac_to_cart_rotation(rotation, lattice):
    """Convert a fractional rotation; lattice vectors are stored as rows."""
    return lattice.T @ rotation @ np.linalg.inv(lattice.T)


def prepare_operations(lattice, positions, rotations, translations, species=None, tolerance=1e-4):
    """Return (0-based atom permutation, Cartesian rotation) for each operation.

    Matching uses periodic Cartesian distances in Angstrom and preserves species.
    The result is geometry-only and can be reused at every photon energy.
    """
    positions = np.asarray(positions, dtype=float)
    rotations = np.asarray(rotations, dtype=float)
    translations = np.asarray(translations, dtype=float)
    natoms = len(positions)
    if (rotations.ndim != 3 or rotations.shape[1:] != (3, 3) or not len(rotations)
            or translations.shape != (len(rotations), 3)):
        raise ValueError('Invalid space-group operations')
    symbols = None if species is None else np.asarray(species)
    if symbols is not None and len(symbols) != natoms:
        raise ValueError('Species length differs from atomic structure')
    operations = []
    for fractional, translation in zip(rotations, translations):
        cartesian = frac_to_cart_rotation(fractional, lattice)
        if not np.allclose(cartesian @ cartesian.T, np.eye(3), atol=1e-8, rtol=0):
            raise ValueError('Space-group Cartesian rotation must be orthogonal')
        permutation = []
        for atom, moved in enumerate(positions @ fractional.T + translation):
            delta = moved - positions
            delta -= np.rint(delta)
            distances = np.linalg.norm(delta @ lattice, axis=1)
            if symbols is not None:
                distances = np.where(symbols == symbols[atom], distances, np.inf)
            target = int(np.argmin(distances))
            if not np.isfinite(distances[target]) or distances[target] > tolerance:
                raise ValueError('Space-group operation does not bijectively map the atomic structure')
            permutation.append(target)
        if len(set(permutation)) != natoms:
            raise ValueError('Space-group operation does not bijectively map the atomic structure')
        operations.append((np.asarray(permutation), cartesian))
    return operations


@dataclass(frozen=True)
class CrystalSymmetry:
    fractional_rotations: np.ndarray
    translations: np.ndarray
    atom_mapping: dict
    sites: dict
    operations: list
    transfers: dict


def _site_entries(data):
    return {int(entry['atom']): {
        'wyckoff': entry.get('Wyckoff'), 'point_group': entry.get('site_point_group'),
        'rotations': np.asarray(entry['rotations'], dtype=float),
    } for entry in data.get('site_symmetries', [])}


def read_symmetry(structure, path='symmetry'):
    """Read once and prepare the symmetry shared by atomic and mode workflows."""
    with open(path) as stream:
        data = yaml.safe_load(stream)
    rotations = np.asarray([op['rotation'] for op in data['space_group_operations']])
    translations = np.asarray([op['translation'] for op in data['space_group_operations']])
    mapping = {int(k): int(v) for k, v in data['atom_mapping'].items()}
    operations = prepare_operations(structure.lattice, structure.fractional_positions,
                                    rotations, translations, structure.atom_symbols)
    chosen = _transfer_indices(mapping, rotations, operations)
    transfers = {pair: operations[i][1] for pair, i in chosen.items()}
    return CrystalSymmetry(rotations, translations, mapping, _site_entries(data), operations, transfers)


def symmetry_operations(structure, symmetry_path='symmetry', tolerance=1e-4):
    """Read operations for mode analysis, including files without site metadata."""
    with open(symmetry_path) as stream:
        data = yaml.safe_load(stream)
    return prepare_operations(structure.lattice, structure.fractional_positions,
                              [op['rotation'] for op in data['space_group_operations']],
                              [op['translation'] for op in data['space_group_operations']],
                              structure.atom_symbols, tolerance)


def read_symmetry_file(filepath="symmetry"):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    rotations = np.array([op['rotation'] for op in data['space_group_operations']])
    translations = np.array([op['translation'] for op in data['space_group_operations']])
    atom_mapping = {int(k): int(v) for k, v in data['atom_mapping'].items()}

    return rotations, translations, atom_mapping


def read_site_symmetries(filepath="symmetry"):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    return _site_entries(data)

def _transfer_indices(atom_mapping, rotations, operations):
    """Choose the same operation for both fractional and Cartesian consumers."""
    chosen = {}
    for atom, representative in atom_mapping.items():
        candidates = [i for i, (permutation, _) in enumerate(operations)
                      if permutation[representative - 1] == atom - 1]
        if not candidates:
            raise ValueError(f'No symmetry operation connects atoms {representative} and {atom}')
        # Preserve the diagonal-first preference and operation ordering.
        chosen[(representative, atom)] = next(
            (i for i in candidates if np.array_equal(rotations[i], np.diag(np.diag(rotations[i])))),
            candidates[0])
    return chosen


def compute_mapping_matrices(positions, rotations, translations, atom_mapping, lattice, species=None):
    """Select fractional transfers using physical-lattice atom permutations."""
    operations = prepare_operations(lattice, positions, rotations, translations, species)
    representatives = sorted(set(atom_mapping.values()))
    groups = {atom: [a for a, rep in atom_mapping.items() if rep == atom] for atom in representatives}
    chosen = _transfer_indices(atom_mapping, rotations, operations)
    return representatives, groups, {pair: rotations[i] for pair, i in chosen.items()}


def run_mapping():
    required_files = ["CONTCAR", "symmetry"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    print("Reading CONTCAR and symmetry files...")
    structure = read_structure("CONTCAR")
    symmetry = read_symmetry(structure)

    print("Finding symmetry matrices that map equivalent atoms...")
    inequivalent_indices = sorted(set(symmetry.atom_mapping.values()))
    equivalent_groups = {rep: [atom for atom, mapped in symmetry.atom_mapping.items() if mapped == rep]
                         for rep in inequivalent_indices}
    inverse = np.linalg.inv(structure.lattice.T)
    final_mapping_matrices = {pair: inverse @ rotation @ structure.lattice.T
                              for pair, rotation in symmetry.transfers.items()}
    num_inequivalent = len(inequivalent_indices)
    print(f"Found {num_inequivalent} inequivalent atoms.")

    print("Writing symmetry_operation_matrices file...")
    with open("symmetry_operation_matrices", "w") as f:
        f.write(f"Number_of_symmetry_independent_atoms:   {num_inequivalent}\n")
        f.write("Indices_of_symmetry_independent_atoms: ")
        f.write(" ".join(map(str, inequivalent_indices)) + "\n")

        for ineq_idx in inequivalent_indices:
            equiv_atoms = sorted(equivalent_groups[ineq_idx])
            f.write(f"Number_of_symmetry_equivalent_atoms_for_atom {ineq_idx} is   {len(equiv_atoms)}\n")
            f.write("Their_indices_are: ")
            f.write(" ".join(map(str, equiv_atoms)) + "\n")
        
        for (ineq_idx, eq_idx), matrix in sorted(final_mapping_matrices.items()):
            f.write(f"\nFind the symmetry operation matrix between atom {ineq_idx} and atom {eq_idx}\n")
            for row in matrix:
                f.write(f"{row[0]:15.4f}{row[1]:15.4f}{row[2]:15.4f}\n")

    print("\nspectropy symmetry finished successfully.")

if __name__ == "__main__":
    run_mapping()
