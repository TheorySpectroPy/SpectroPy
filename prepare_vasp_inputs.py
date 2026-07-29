import numpy as np
import os
import sys
import shutil
from collections import defaultdict

def run_generate_displacements():
    """
    Reads CONTCAR and displacements.dat to generate displaced POSCAR directories.
    """
    # --- Check for required input files ---
    if not os.path.exists("CONTCAR"):
        print("***** CONTCAR not found *****")
        sys.exit(1)
    if not os.path.exists("displacements.dat"):
        print("***** displacements.dat not found *****")
        sys.exit(1)
        
    print("Generates all POSCARs for VASP needed to find dielectric tensors.")
    print("Necessary displacements must be given in displacements.dat.")
    print("This program reads CONTCAR.")
    
    print("-" * 40)

    # --- 1. Read CONTCAR to get original structure data ---
    print("Reading CONTCAR...")
    with open("CONTCAR", 'r') as f:
        lines = f.readlines()

    comment = lines[0]
    lattice_constant_line = lines[1]
    lattice_vectors_lines = lines[2:5]

    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        atom_symbols_line = lines[5]
        coord_type_line_idx = 7
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        # Create a placeholder for atom symbols line
        atom_symbols_line = " ".join([f"El{i+1}" for i in range(len(atom_counts))]) + "\n"
        coord_type_line_idx = 6

    total_atoms = np.sum(atom_counts)
    coord_type_line = lines[coord_type_line_idx]
    
    # Store original positions in fractional coordinates
    original_positions = np.array([
        list(map(float, line.split()[:3])) 
        for line in lines[coord_type_line_idx+1 : coord_type_line_idx+1+total_atoms]
    ])

    # --- 2. Copy CONTCAR to ref_poscar.vasp ---
    print("Copying CONTCAR to ref_poscar.vasp...")
    shutil.copyfile("CONTCAR", "ref_poscar.vasp")

    # --- 3. Read displacements.dat to get displacement data ---
    print("Reading displacements.dat...")
    with open("displacements.dat", 'r') as f:
        displacements_lines = f.readlines()
    
    num_disps = int(displacements_lines[1].split()[0])
    
    displacements = []
    for line in displacements_lines[2 : 2 + num_disps]:
        parts = line.split()
        label = parts[0]
        atom_index = int(parts[1]) # This is 1-based
        vector = np.array(list(map(float, parts[2:5])))
        displacements.append({'label': label, 'index': atom_index, 'vector': vector})

    # --- 4. Generate displaced POSCAR files and calculation directories ---
    print(f"Generating {num_disps} displaced POSCAR files (poscar*)...")
    
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    # Use a dictionary to count displacements for each atom label (e.g., W1, Te1)
    disp_counter = defaultdict(int)

    for disp_data in displacements:
        label = disp_data['label']
        atom_idx_1based = disp_data['index']
        disp_vector = disp_data['vector']

        # Generate the unique suffix (e.g., 'a', 'b', 'c'...).
        suffix = alphabet[disp_counter[label]]
        disp_counter[label] += 1
        poscar_filename = f"pos_atom{atom_idx_1based}{suffix}"
        calculation_dir = f"ra_{poscar_filename}"

        # Write a standalone POSCAR and put an identical copy in the DFT
        # directory. No INCAR/KPOINTS/POTCAR links or shell scripts are made:
        # the user supplies their own DFT inputs and execution method.
        with open(poscar_filename, "w") as f_pos:
            f_pos.write(comment)
            f_pos.write(lattice_constant_line)
            f_pos.writelines(lattice_vectors_lines)
            f_pos.write(atom_symbols_line)
            f_pos.write(" ".join(map(str, atom_counts)) + "\n")
            f_pos.write(coord_type_line)

            new_positions = np.copy(original_positions)
            new_positions[atom_idx_1based - 1] += disp_vector
            for pos in new_positions:
                f_pos.write(f"  {pos[0]:.16f} {pos[1]:.16f} {pos[2]:.16f}\n")

        os.makedirs(calculation_dir, exist_ok=True)
        shutil.copyfile(poscar_filename, os.path.join(calculation_dir, "POSCAR"))
        
    print("\nprepare_vasp_inputs.py finished successfully.")
    print("Generated 'pos_atom*' files and all 'ra_pos_atom*' directories.")
    print("\nAdd your DFT inputs and run your calculations in each directory.")

# This makes the script runnable from the command line
if __name__ == "__main__":
    run_generate_displacements()
