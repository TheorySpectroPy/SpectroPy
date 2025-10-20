import numpy as np
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

def read_static_diel_from_xml(xml_path):
    """
    Parses a vasprun.xml file to find the static dielectric tensor.

    Args:
        xml_path (str): Path to the vasprun.xml file.

    Returns:
        np.ndarray or None: The 3x3 real dielectric tensor, or None if parsing fails.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except (ET.ParseError, FileNotFoundError):
        print(f"Error: Could not parse or find {xml_path}")
        return None

    try:
        # VASP writes the macroscopic static tensor under the "epsilon" or "epsilon_scf" tag
        diel_tensor_element = root.find("./calculation/dielectric/varray[@name='epsilon']")
        if diel_tensor_element is None:
            diel_tensor_element = root.find("./calculation/dielectric/varray[@name='epsilon_scf']")
        
        if diel_tensor_element is None:
            raise AttributeError # Trigger the except block if not found
            
        rows = []
        for r in diel_tensor_element.findall('r'):
            rows.append(list(map(float, r.text.split())))
        
        return np.array(rows)

    except AttributeError:
        print(f"Error: Could not find static dielectric tensor in {xml_path}")
        return None

def run_calculate_dielectric_derivatives_static():

    # --- Check for required input files ---
    required_files = ["ref_poscar.vasp", "dielectric_derivatives.dat"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)
    if not os.path.exists("vasprun"):
        print("***** vasprun directory not found. Did you run 'kopia'? *****")
        sys.exit(1)
        
    # --- Preamble ---
    print("Program <calculate_dielectric_derivatives_static.py>")
    print("Calculates STATIC dielectric tensor derivatives from vasprun.xml files.")
    print("-" * 40)

    # --- 1. Read ref_poscar.vasp ---
    print("Reading ref_poscar.vasp...")
    with open("ref_poscar.vasp", 'r') as f:
        lines = f.readlines()
    ref_poscar.vasp_comment = lines[0].strip()
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    
    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[8:8+total_atoms]])
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[7:7+total_atoms]])

    # --- 2. Read dielectric_derivatives.dat ---
    print("Reading dielectric_derivatives.dat...")
    with open("dielectric_derivatives.dat", 'r') as f:
        displacements_lines = f.readlines()
    num_disps = int(displacements_lines[1].split()[0])
    
    displacements = []
    for line in displacements_lines[2 : 2 + num_disps]:
        parts = line.split()
        displacements.append({'label': parts[0], 'index': int(parts[1]), 'vector': np.array(list(map(float, parts[2:5])))})
    
    # --- 3. Process XML files and collect dielectric tensors ---
    print("Processing vasprun.xml files from vasprun directory...")
    all_epsilons = np.zeros((num_disps, 3, 3))
    cart_disps = np.zeros((num_disps, 3))
    
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    disp_counter = defaultdict(int)

    # Open log file
    with open("RAMAN.out", 'w') as f_log:
        f_log.write("Dielectric constant tensors from the following files:\n")
        f_log.write("Atom position in fractional coordinates with respect to matrix:\n")
        for vec in lattice_vectors:
            f_log.write(f"   {vec[0]:21.16f} {vec[1]:21.16f} {vec[2]:21.16f}\n")
        f_log.write("-" * 80 + "\n")
        
        for i, disp_data in enumerate(displacements):
            label = disp_data['label']
            suffix = alphabet[disp_counter[label]]
            disp_counter[label] += 1
            
            xml_filename = f"{label}{suffix}.xml"
            xml_path = os.path.join("vasprun", xml_filename)
            
            epsilon = read_static_diel_from_xml(xml_path)
            if epsilon is None:
                sys.exit(1) # Error message is printed inside the function
                
            all_epsilons[i] = epsilon
            
            # Convert fractional displacement vector to Cartesian
            cart_disp_vec = disp_data['vector'] @ lattice_vectors
            cart_disps[i] = cart_disp_vec
            
            # Write to log file
            pos_str = " ".join([f"{p:8.5f}" for p in positions[disp_data['index']-1]])
            disp_str = " ".join([f"{d:8.5f}" for d in cart_disp_vec])
            f_log.write(f"{xml_filename:<12s}   {pos_str}      {disp_str}\n")
            for j in range(3):
                row_str = " ".join([f"{x:12.6f}" for x in epsilon[j]])
                f_log.write(f"  {row_str}\n")
            f_log.write("\n")

    print("Finished processing XMLs. Log written to RAMAN.out")

    # --- 4. Calculate Raman Tensors via Finite Difference ---
    print("Calculating Raman tensors using finite difference method...")
    num_ineq_atoms = num_disps // 6
    # d(eps)/du_alpha for each inequivalent atom
    eps_deriv = np.zeros((num_ineq_atoms, 3, 3, 3)) # (atom, alpha, i, j)

    for i in range(num_ineq_atoms):
        idx_base = i * 6
        for alpha in range(3):
            idx_plus = idx_base + alpha * 2
            idx_minus = idx_base + alpha * 2 + 1
            
            eps_plus = all_epsilons[idx_plus]
            eps_minus = all_epsilons[idx_minus]
            
            denominator = np.linalg.norm(cart_disps[idx_plus] - cart_disps[idx_minus])
            if denominator < 1e-6: continue
            
            eps_deriv[i, alpha, :, :] = (eps_plus - eps_minus) / denominator
    
    
    pi = np.pi
    cell_volume = np.linalg.det(lattice_vectors)
    conv_factor = cell_volume / (4 * pi)
    eps_deriv *= conv_factor

    with open("dielectric_derivatives.dat", 'w') as f:
        f.write(f"! Static Raman Tensors calculated from vasprun.xml files\n")
        f.write("!\n! Unit cell matrix:\n")
        for vec in lattice_vectors:
            f.write(f"!   {vec[0]:21.16f} {vec[1]:21.16f} {vec[2]:21.16f}\n")
        f.write("!\n! DIELECTRIC DERIVATIVES\n")
        
        unique_atom_labels = [d['label'] for d in displacements[::6]]
        unique_atom_indices = [d['index'] for d in displacements[::6]]

        for i, label in enumerate(unique_atom_labels):
            pos = positions[unique_atom_indices[i]-1]
            f.write(f"      {label:<4s} {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}\n")
            # This is a 3x9 tensor: (d(eps)/dx, d(eps)/dy, d(eps)/dz)
            # Reshape for printing:
            tensor_to_print = eps_deriv[i].transpose(0, 2, 1).reshape(3, 9)
            for row in tensor_to_print:
                f.write("".join([f"{x:16.4f}" for x in row]) + "\n")

    print("\ncalculate_dielectric_derivatives_static.py finished successfully.")

if __name__ == "__main__":
    run_calculate_dielectric_derivatives_static()