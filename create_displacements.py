import numpy as np
import os
import sys

def run_displacements(contcar_path="CONTCAR"):
    """
    Reads a VASP CONTCAR file and generates the necessary displacement files
    for a Raman calculation.

    Args:
        contcar_path (str): The path to the CONTCAR file.
    """
    # Check if CONTCAR exists
    if not os.path.exists(contcar_path):
        print(f"Error: {contcar_path} not found.")
        sys.exit(1)

    # --- 1. Read CONTCAR ---
    print(f"Reading structure from {contcar_path}...")
    with open(contcar_path, 'r') as f:
        lines = f.readlines()

    comment = lines[0].strip()
    lattice_constant = float(lines[1].strip())
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    lattice_vectors *= lattice_constant

    # Determine if atom symbols are present (VASP 5 format)
    try:
        # Check if the 6th line contains numbers (atom counts)
        atom_counts = np.array(lines[6].split(), dtype=int)
        atom_symbols = lines[5].split()
        coord_type_line = 7
    except (ValueError, IndexError):
        # If not, assume VASP 4 format (no symbols line)
        atom_counts = np.array(lines[5].split(), dtype=int)
        atom_symbols = ["Elem" + str(i+1) for i in range(len(atom_counts))] # Placeholder symbols
        coord_type_line = 6
        print("Warning: Atom symbols not found in CONTCAR, using placeholders.")

    total_atoms = np.sum(atom_counts)
    coord_type = lines[coord_type_line].strip()
    
    # Read positions
    positions_raw = np.array([list(map(float, line.split()[:3])) for line in lines[coord_type_line+1 : coord_type_line+1+total_atoms]])

    # Convert to Cartesian if necessary
    if coord_type.lower().startswith('d'): # Direct coordinates
        positions_cart = positions_raw @ lattice_vectors
        positions_frac = positions_raw
    else: # Cartesian coordinates
        positions_cart = positions_raw
        positions_frac = positions_cart @ np.linalg.inv(lattice_vectors)

    # Generate a full list of atom symbols
    full_atom_symbols = []
    for symbol, count in zip(atom_symbols, atom_counts):
        full_atom_symbols.extend([symbol] * count)

    # --- 2. Define Cartesian Displacements ---
    displacement_val = 0.01  # Hardcoded value for now
    cart_disps = np.array([
        [ displacement_val,  0.0,  0.0],
        [-displacement_val,  0.0,  0.0],
        [ 0.0,  displacement_val,  0.0],
        [ 0.0, -displacement_val,  0.0],
        [ 0.0,  0.0,  displacement_val],
        [ 0.0,  0.0, -displacement_val]
    ])

    # --- 3. Convert Displacements to Fractional Coords ---
    print("Calculating fractional displacements...")
    inv_lattice = np.linalg.inv(lattice_vectors)
    frac_disps = cart_disps @ inv_lattice

    # --- 4. Write ref_poscar.vasp (equilibrium structure) ---
    print("Writing ref_poscar.vasp...")
    with open("ref_poscar.vasp", "w") as f:
        f.write(f"{comment}\n")
        f.write("   1.0\n")
        for vec in lattice_vectors:
            f.write(f"   {vec[0]:21.16f} {vec[1]:21.16f} {vec[2]:21.16f}\n")
        f.write(f"   {' '.join(atom_symbols)}\n")
        f.write(f"   {' '.join(map(str, atom_counts))}\n")
        f.write("Direct\n")
        for pos in positions_frac:
            f.write(f"   {pos[0]:19.16f} {pos[1]:19.16f} {pos[2]:19.16f}\n")

    # --- 5. Write displacements.dat (displacement list) ---
    print("Writing displacements.dat...")
    # Create unique atom labels like 'Si1', 'Si2', etc.
    atom_counts_dict = {symbol: 0 for symbol in atom_symbols}
    ineq_atom_symbols = []
    for symbol in full_atom_symbols:
        atom_counts_dict[symbol] += 1
        ineq_atom_symbols.append(f"{symbol}{atom_counts_dict[symbol]}")

    with open("displacements.dat", "w") as f:
        f.write("displacements for VASP. System\n")
        f.write(f"{total_atoms * 6:6d}   Number of displacements\n")
        
        # 1 indexed
        for atom_idx in range(total_atoms):
            atom_label = ineq_atom_symbols[atom_idx]
            for disp in frac_disps:
                f.write(f"{atom_label:<5s}{atom_idx + 1:6d}{disp[0]:12.6f}{disp[1]:12.6f}{disp[2]:12.6f}\n")
        
        f.write(f"{total_atoms:6d}{total_atoms:6d}     Number of atoms in SC\n")
        for atom_idx in range(total_atoms):
            atom_label = ineq_atom_symbols[atom_idx]
            f.write(f"{atom_label:<5s}{atom_idx + 1:6d}\n")

    print("\ndisplacements.dat.py finished successfully.")
    print("Generated ref_poscar.vasp and displacements.dat.")

if __name__ == "__main__":
    run_displacements()