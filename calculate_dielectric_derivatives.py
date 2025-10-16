import numpy as np
import os
import sys
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict

def read_diel_from_xml(xml_path, target_frequency):
    """
    Parses a vasprun.xml file to find the complex dielectric tensor
    at the frequency closest to the target frequency.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except (ET.ParseError, FileNotFoundError):
        print(f"Error: Could not parse or find {xml_path}")
        return None, None

    try:
        diel_func = root.find("./calculation/dielectricfunction")
        imag_part = diel_func.find('./imag/array/set')
        real_part = diel_func.find('./real/array/set')

        frequencies = []
        eps_imag_components = []
        eps_real_components = []

        for r in imag_part.findall('r'):
            data = list(map(float, r.text.split()))
            frequencies.append(data[0])
            eps_imag_components.append(data[1:])
        
        for r in real_part.findall('r'):
            data = list(map(float, r.text.split()))
            eps_real_components.append(data[1:])

        frequencies = np.array(frequencies)
        eps_imag_components = np.array(eps_imag_components)
        eps_real_components = np.array(eps_real_components)
        
        idx = np.argmin(np.abs(frequencies - target_frequency))
        
        ri = eps_real_components[idx]
        im = eps_imag_components[idx]
        
        eps_real = np.array([[ri[0], ri[3], ri[5]], [ri[3], ri[1], ri[4]], [ri[5], ri[4], ri[2]]])
        eps_imag = np.array([[im[0], im[3], im[5]], [im[3], im[1], im[4]], [im[5], im[4], im[2]]])

        return eps_real, eps_imag

    except AttributeError:
        print(f"Error: Could not find dielectric data in {xml_path}")
        return None, None

def run_generate_derivatives():
    """
    Main function to collect VASP outputs and calculate derivatives of dielectric tensors.
    """
    required_files = ["ref_poscar.vasp", "displacements.dat"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)
    
    print("Reading displacements.dat to determine which files to collect...")
    with open("displacements.dat", 'r') as f:
        displacements_lines = f.readlines()
    num_disps = int(displacements_lines[1].split()[0])
    
    displacements = [{'label': line.split()[0]} for line in displacements_lines[2 : 2 + num_disps]]

    print("Collecting vasprun.xml files...")
    os.makedirs("AXML", exist_ok=True)
    
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    disp_counter = defaultdict(int)
    
    for disp_data in displacements:
        label = disp_data['label']
        suffix = alphabet[disp_counter[label]]
        disp_counter[label] += 1
        
        source_dir_name = f"ra_pos_{label}{suffix}"
        source_path = os.path.join(source_dir_name, "vasprun.xml")
        dest_filename = f"{label}{suffix}.xml"
        dest_path = os.path.join("AXML", dest_filename)
        
        try:
            shutil.copyfile(source_path, dest_path)
        except FileNotFoundError:
            print(f"Error: Source file not found: {source_path}")
            print("Please ensure all VASP calculations have finished.")
            sys.exit(1)
            
    print("File collection complete. AXML directory is ready.")
    print("-" * 40)
        
    print("Program <calculate_dielectric_derivatives.py>")
    print("Calculates derivatives of dielectric tensors from the collected vasprun.xml files.")
    try:
        target_frequency = float(input("Please input the laser frequency (eV): "))
    except ValueError:
        print("Invalid frequency. Exiting.")
        sys.exit(1)
    print("-" * 40)

    print("Reading ref_poscar.vasp...")
    with open("ref_poscar.vasp", 'r') as f:
        lines = f.readlines()
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    
    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[8:8+total_atoms]])
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[7:7+total_atoms]])

    displacements = []
    for line in displacements_lines[2 : 2 + num_disps]:
        parts = line.split()
        displacements.append({'label': parts[0], 'index': int(parts[1]), 'vector': np.array(list(map(float, parts[2:5])))})
    
    print("Processing vasprun.xml files from AXML directory...")
    all_eps_real = np.zeros((num_disps, 3, 3))
    all_eps_imag = np.zeros((num_disps, 3, 3))
    cart_disps = np.zeros((num_disps, 3))
    
    disp_counter = defaultdict(int)

    epsilon_log_filename = f"EPSILON_{target_frequency:.2f}.dat"
    with open(epsilon_log_filename, 'w') as f_log:
        f_log.write(f"Dielectric tensors for target frequency {target_frequency:.4f} eV\n")
        f_log.write("File                 Positions (frac)              Displacement (Angstrom)\n")
        f_log.write("-" * 80 + "\n")
        
        for i, disp_data in enumerate(displacements):
            label = disp_data['label']
            suffix = alphabet[disp_counter[label]]
            disp_counter[label] += 1
            
            xml_filename = f"{label}{suffix}.xml"
            xml_path = os.path.join("AXML", xml_filename)
            
            eps_real, eps_imag = read_diel_from_xml(xml_path, target_frequency)
            if eps_real is None: sys.exit(1)
            
            all_eps_real[i] = eps_real
            all_eps_imag[i] = eps_imag
            cart_disp_vec = disp_data['vector'] @ lattice_vectors
            cart_disps[i] = cart_disp_vec
            
            pos_str = " ".join([f"{p:8.5f}" for p in positions[disp_data['index']-1]])
            disp_str = " ".join([f"{d:8.5f}" for d in cart_disp_vec])
            f_log.write(f"{xml_filename:<12s}   {pos_str}      {disp_str}\n")
            for j in range(3):
                real_str = " ".join([f"{x:10.6f}" for x in eps_real[j]])
                imag_str = " ".join([f"{x:10.6f}" for x in eps_imag[j]])
                f_log.write(f"  real: {real_str} | imag: {imag_str}\n")
            f_log.write("\n")

    print(f"Finished processing XMLs. Log written to {epsilon_log_filename}")

    print("Calculating derivatives of dielectric tensors using finite difference method...")
    num_ineq_atoms = num_disps // 6
    eps_deriv_real = np.zeros((num_ineq_atoms, 3, 3, 3))
    eps_deriv_imag = np.zeros((num_ineq_atoms, 3, 3, 3))

    for i in range(num_ineq_atoms):
        idx_base = i * 6
        for alpha in range(3):
            idx_plus = idx_base + alpha * 2
            idx_minus = idx_base + alpha * 2 + 1
            
            denominator = np.linalg.norm(cart_disps[idx_plus] - cart_disps[idx_minus])
            if denominator < 1e-6: continue
            
            eps_deriv_real[i, alpha, :, :] = (all_eps_real[idx_plus] - all_eps_real[idx_minus]) / denominator
            eps_deriv_imag[i, alpha, :, :] = (all_eps_imag[idx_plus] - all_eps_imag[idx_minus]) / denominator
    
    dielectric_derivatives_filename = f"dielectric_derivatives_{target_frequency:.2f}"
    print(f"Writing final derivatives of dielectric tensors to {dielectric_derivatives_filename}...")
    
    pi = np.pi
    cell_volume = np.linalg.det(lattice_vectors)
    conv_factor = cell_volume / (4 * pi)
    
    eps_deriv_real *= conv_factor
    eps_deriv_imag *= conv_factor

    with open(dielectric_derivatives_filename, 'w') as f:
        f.write(f"! derivatives of dielectric tensors calculated for laser frequency {target_frequency:.4f} eV\n")
        f.write("! from VASP vasprun.xml files\n")
        f.write("!\n! Unit cell matrix:\n")

        for vec in lattice_vectors.T:
            f.write(f"!   {vec[0]:21.16f} {vec[1]:21.16f} {vec[2]:21.16f}\n")
        f.write("!\n! derivatives of dielectric tensors\n")
        
        unique_atom_labels = [d['label'] for d in displacements[::6]]
        unique_atom_indices = [d['index'] for d in displacements[::6]]

        f.write("! The Real Part of Raman tensor:\n")
        for i, label in enumerate(unique_atom_labels):
            pos = positions[unique_atom_indices[i]-1]
            f.write(f"      {label:<4s} {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}\n")
            
            tensor_x = eps_deriv_real[i, 0, :, :]
            tensor_y = eps_deriv_real[i, 1, :, :]
            tensor_z = eps_deriv_real[i, 2, :, :]
            tensor_to_print = np.hstack((tensor_x, tensor_y, tensor_z))
            
            for row in tensor_to_print:
                f.write("".join([f"{x:16.4f}" for x in row]) + "\n")
        
        f.write("\n! The Imaginary Part of Raman tensor:\n")
        for i, label in enumerate(unique_atom_labels):
            pos = positions[unique_atom_indices[i]-1]
            f.write(f"      {label:<4s} {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}\n")

            tensor_x_imag = eps_deriv_imag[i, 0, :, :]
            tensor_y_imag = eps_deriv_imag[i, 1, :, :]
            tensor_z_imag = eps_deriv_imag[i, 2, :, :]
            tensor_to_print_imag = np.hstack((tensor_x_imag, tensor_y_imag, tensor_z_imag))

            for row in tensor_to_print_imag:
                f.write("".join([f"{x:16.4f}" for x in row]) + "\n")
                
    print("\ncalculate_dielectric_derivatives.py finished successfully.")

if __name__ == "__main__":
    run_generate_derivatives()