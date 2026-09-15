import numpy as np
import os
import sys
import yaml
from spectropy_structure import read_structure


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

    site_symmetries = {}
    for entry in data['site_symmetries']:
        site_symmetries[int(entry['atom'])] = {
            "wyckoff": entry['Wyckoff'],
            "point_group": entry['site_point_group'],
            "rotations": np.array(entry['rotations']),
        }
    return site_symmetries

def compute_mapping_matrices(positions, rotations, translations, atom_mapping):
    inequivalent_indices = sorted(list(set(atom_mapping.values())))

    equivalent_groups = {idx: [] for idx in inequivalent_indices}
    for atom_idx, maps_to_idx in atom_mapping.items():
        equivalent_groups[maps_to_idx].append(atom_idx)

    final_mapping_matrices = {}

    #Determine which symmetry operation connects the equivalent atom to the inequivalent one
    for ineq_idx in inequivalent_indices:
        for eq_idx in equivalent_groups[ineq_idx]:
            pos_ineq = positions[ineq_idx - 1]
            pos_eq = positions[eq_idx - 1]

            found_matrices = []
            for i in range(len(rotations)):
                rot = rotations[i]
                trans = translations[i]

                pos_prime = rot @ pos_ineq + trans

                #Subtract nearest lattice translation
                delta = pos_prime - pos_eq
                periodic_delta = delta - np.round(delta)

                if np.allclose(periodic_delta, 0, atol=1e-5):
                    found_matrices.append(rot)

            if not found_matrices:
                print(f"Warning: No symmetry operation found between atom {ineq_idx} and {eq_idx}!")
                continue

            #Preferentially pick a diagonal matrix as the representative operation
            chosen_matrix = None
            for matrix in found_matrices:
                if np.count_nonzero(matrix - np.diag(np.diagonal(matrix))) == 0:
                    chosen_matrix = matrix
                    break

            if chosen_matrix is None:
                chosen_matrix = found_matrices[0]

            final_mapping_matrices[(ineq_idx, eq_idx)] = chosen_matrix

    return inequivalent_indices, equivalent_groups, final_mapping_matrices


def run_mapping():
    required_files = ["CONTCAR", "symmetry"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    print("Reading CONTCAR and symmetry files...")
    positions = read_structure("CONTCAR").fractional_positions
    rotations, translations, atom_mapping = read_symmetry_file()

    print("Finding symmetry matrices that map equivalent atoms...")
    inequivalent_indices, equivalent_groups, final_mapping_matrices = compute_mapping_matrices(
        positions, rotations, translations, atom_mapping
    )
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

    print("\nprocess_symmetry.py finished successfully.")

if __name__ == "__main__":
    run_mapping()
