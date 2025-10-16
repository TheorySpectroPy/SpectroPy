import numpy as np
import os
import sys
import shutil
from collections import defaultdict

def run_generate_displacements():
    """
    Reads CONTCAR and displacements.dat to generate a series of displaced POSCAR
    files and a run script (run_vasp_calcs.sh) to manage VASP calculations.
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
    print("It also generates <run_vasp_calcs.sh> script to set up calculations.")
    
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

    # --- 4. Generate displaced POSCAR files and 'collect_results.sh' script ---
    print(f"Generating {num_disps} displaced POSCAR files (pos_*)...")
    
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    # Use a dictionary to count displacements for each atom label (e.g., W1, Te1)
    disp_counter = defaultdict(int)

    with open("collect_results.sh", "w") as f_collect:
        f_collect.write("#!/bin/bash\n")
        f_collect.write("#\n")
        f_collect.write("mkdir -p AXML\n")
        f_collect.write("#\n")

        for disp_data in displacements:
            label = disp_data['label']
            atom_idx_1based = disp_data['index']
            disp_vector = disp_data['vector']

            # Generate the unique suffix (e.g., 'a', 'b', 'c'...)
            suffix = alphabet[disp_counter[label]]
            disp_counter[label] += 1
            
            poscar_filename = f"pos_{label}{suffix}"
            
            # Write the copy command to the collect_results.sh script
            xml_name = f"{label}{suffix}"
            f_collect.write(f"cp ra_{poscar_filename}/vasprun.xml AXML/{xml_name}.xml\n")

            # Create the new displaced POSCAR
            with open(poscar_filename, "w") as f_pos:
                f_pos.write(comment)
                f_pos.write(lattice_constant_line)
                f_pos.writelines(lattice_vectors_lines)
                f_pos.write(atom_symbols_line)
                f_pos.write(" ".join(map(str, atom_counts)) + "\n")
                f_pos.write(coord_type_line)

                # Create a copy of positions to modify
                new_positions = np.copy(original_positions)
                
                # Apply displacement (convert 1-based index to 0-based)
                atom_idx_0based = atom_idx_1based - 1
                new_positions[atom_idx_0based] += disp_vector
                
                for pos in new_positions:
                    f_pos.write(f"  {pos[0]:.16f} {pos[1]:.16f} {pos[2]:.16f}\n")

    # --- 5. Generate the 'setup_vasp_calcs.sh' script ---
    print("Generating setup_vasp_calcs.sh script...")
    with open("setup_vasp_calcs.sh", "w") as f_run:
        f_run.write("#!/bin/bash\n")
        f_run.write("# This script sets up a calculation directory for each displaced POSCAR.\n")
        f_run.write("# After running this, you will need to run VASP in each 'ra_pos_*' subdirectory.\n")
        f_run.write("\n")
        f_run.write("for d in pos_*; do\n")
        f_run.write("  echo \"Setting up directory for $d\"\n")
        f_run.write("  mkdir -p ra_$d\n")
        f_run.write("  cp $d ra_$d/POSCAR\n")
        f_run.write("  ln -sf ../KPOINTS ra_$d/KPOINTS\n")
        f_run.write("  ln -sf ../POTCAR ra_$d/POTCAR\n")
        f_run.write("  ln -sf ../INCAR ra_$d/INCAR\n")
        if os.path.exists("vdw_kernel.bindat"):
            f_run.write("  ln -sf ../vdw_kernel.bindat ra_$d/vdw_kernel.bindat\n")
        f_run.write("done\n")
        f_run.write("\n")
        f_run.write("echo \"All calculation directories have been set up.\"\n")
        
    print("\nprepare_vasp_inputs.py finished successfully.")
    print("Generated 'pos_*' files, 'collect_results.sh', and 'setup_vasp_calcs.sh'.")
    print("Next steps:")
    print("1. Make the setup script executable: chmod +x setup_vasp_calcs.sh")
    print("2. Execute the script to create directories: ./setup_vasp_calcs.sh")

# This makes the script runnable from the command line
if __name__ == "__main__":
    run_generate_displacements()