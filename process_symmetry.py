import numpy as np
import os
import sys
import yaml

def read_contcar(filepath="CONTCAR"):
    """Reads a VASP CONTCAR file to get atomic positions."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    lattice_constant = float(lines[1].strip())
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    lattice_vectors *= lattice_constant

    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        coord_type_line_idx = 7
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        coord_type_line_idx = 6

    total_atoms = np.sum(atom_counts)
    coord_type = lines[coord_type_line_idx].strip()
    
    positions_raw = np.array([list(map(float, line.split()[:3])) for line in lines[coord_type_line_idx+1 : coord_type_line_idx+1+total_atoms]])

    if coord_type.lower().startswith('d'):
        return positions_raw, total_atoms
    else: # Convert Cartesian to Direct/Fractional
        inv_lattice = np.linalg.inv(lattice_vectors)
        return positions_raw @ inv_lattice, total_atoms

def read_symmetry_file(filepath="symmetry"):
    """
    Parses phonopy's `--symmetry` output file (real YAML) to extract the
    space-group operations and the atom-equivalence mapping.

    The previous version of this function hand-parsed the file line by line
    with `next(f)` calls after the `with open(...) as f:` block had already
    closed the handle -- a `ValueError: I/O operation on closed file` on
    every call, so nothing downstream of it ever actually ran. The file is
    valid YAML (phonopy writes it that way), so just load it properly.
    """
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    rotations = np.array([op['rotation'] for op in data['space_group_operations']])
    translations = np.array([op['translation'] for op in data['space_group_operations']])
    atom_mapping = {int(k): int(v) for k, v in data['atom_mapping'].items()}

    return rotations, translations, atom_mapping


def read_site_symmetries(filepath="symmetry"):
    """
    Parses the per-atom `site_symmetries` block phonopy also writes into the
    same file: for each symmetry-inequivalent atom, the subgroup of the space
    group that maps that atom to itself (in fractional coordinates). This is
    what's needed to find the minimal set of displacement directions for
    each atom (see raman_symmetry_utils.find_minimal_directions) -- it is
    *not* the same information as `atom_mapping`/`read_symmetry_file`, which
    only relates *different* atoms to each other, not an atom to itself.

    Returns: {atom_index (1-based): {"wyckoff": str, "point_group": str,
                                      "rotations": (n,3,3) int array}}
    """
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
    """
    Core of run_mapping(), factored out so other scripts (e.g.
    reconstruct_dielectric_derivatives.py's atom-expansion step) can get the
    {(rep_idx, atom_idx): R_frac} dict directly instead of round-tripping
    through the symmetry_operation_matrices text file.

    Returns (inequivalent_indices, equivalent_groups, final_mapping_matrices)
    with final_mapping_matrices keyed (rep_idx, atom_idx) both 1-based,
    R_frac in fractional coordinates (convert with frac_to_cart_rotation
    from reconstruct_dielectric_derivatives.py before applying to a
    Cartesian tensor).
    """
    inequivalent_indices = sorted(list(set(atom_mapping.values())))

    equivalent_groups = {idx: [] for idx in inequivalent_indices}
    for atom_idx, maps_to_idx in atom_mapping.items():
        equivalent_groups[maps_to_idx].append(atom_idx)

    final_mapping_matrices = {}

    for ineq_idx in inequivalent_indices:
        for eq_idx in equivalent_groups[ineq_idx]:
            pos_ineq = positions[ineq_idx - 1] # Convert to 0-based index
            pos_eq = positions[eq_idx - 1]     # Convert to 0-based index

            found_matrices = []
            for i in range(len(rotations)):
                rot = rotations[i]
                trans = translations[i]

                # Apply symmetry operation: r' = R*r + t
                pos_prime = rot @ pos_ineq + trans

                # Check if r' is equivalent to r_eq under periodic boundary conditions
                delta = pos_prime - pos_eq
                periodic_delta = delta - np.round(delta) # Brings components to [-0.5, 0.5]

                if np.allclose(periodic_delta, 0, atol=1e-5):
                    found_matrices.append(rot)

            if not found_matrices:
                print(f"Warning: No symmetry operation found between atom {ineq_idx} and {eq_idx}!")
                continue

            # prefer a diagonal matrix if multiple are found
            chosen_matrix = None
            for matrix in found_matrices:
                # Check if off-diagonal elements are all zero
                if np.count_nonzero(matrix - np.diag(np.diagonal(matrix))) == 0:
                    chosen_matrix = matrix
                    break

            # If no diagonal matrix was found, just pick the first one
            if chosen_matrix is None:
                chosen_matrix = found_matrices[0]

            final_mapping_matrices[(ineq_idx, eq_idx)] = chosen_matrix

    return inequivalent_indices, equivalent_groups, final_mapping_matrices


def run_mapping():
    """Main function to find and write symmetry operation matrices."""
    # --- Check for required input files ---
    required_files = ["CONTCAR", "symmetry"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    # --- Read Input Files ---
    print("Reading CONTCAR and symmetry files...")
    positions, total_atoms = read_contcar()
    rotations, translations, atom_mapping = read_symmetry_file()

    print("Finding symmetry matrices that map equivalent atoms...")
    inequivalent_indices, equivalent_groups, final_mapping_matrices = compute_mapping_matrices(
        positions, rotations, translations, atom_mapping
    )
    num_inequivalent = len(inequivalent_indices)
    print(f"Found {num_inequivalent} inequivalent atoms.")

    # --- Write the Output File ---
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