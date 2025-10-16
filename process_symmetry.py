import numpy as np
import os
import sys

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
    Parses a phonopy-style symmetry file to extract operations and atom mappings.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    rotations = []
    translations = []
    atom_mapping = {}
    
    mode = None
    for line in lines:
        if "space_group_operations" in line:
            mode = "operations"
            continue
        elif "atom_mapping" in line:
            mode = "mapping"
            continue
            
        if mode == "operations" and "rotation" in line:
            rotation_matrix = []
            for _ in range(3):

                line = next(f)
                row = [int(line[3:6]), int(line[6:10]), int(line[10:14])]
                rotation_matrix.append(row)
            rotations.append(np.array(rotation_matrix))
            
            line = next(f) # translation line
            trans_vec = np.array(list(map(float, line.split(":")[1].split())))
            translations.append(trans_vec)

        if mode == "mapping" and ":" in line:
            parts = line.split(":")
            atom_idx = int(parts[0].strip())
            mapped_to_idx = int(parts[1].strip())
            atom_mapping[atom_idx] = mapped_to_idx

    return np.array(rotations), np.array(translations), atom_mapping

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
    
    # --- Process Atom Mappings ---
    # Find the set of unique (inequivalent) atom indices
    inequivalent_indices = sorted(list(set(atom_mapping.values())))
    num_inequivalent = len(inequivalent_indices)

    # Group equivalent atoms under their inequivalent parent
    equivalent_groups = {idx: [] for idx in inequivalent_indices}
    for atom_idx, maps_to_idx in atom_mapping.items():
        equivalent_groups[maps_to_idx].append(atom_idx)

    print(f"Found {num_inequivalent} inequivalent atoms.")
    
    # --- Find the Mapping Matrix for Each Pair ---
    print("Finding symmetry matrices that map equivalent atoms...")
    # This will store the final chosen rotation matrix for each pair
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