import numpy as np
import os
import sys
import yaml
import re

def get_user_input():
    """Gets polarization vectors from 'input' file or user prompt."""
    if os.path.exists("input"):
        print("Reading experimental geometry from 'input' file...")
        settings = np.loadtxt("input", usecols=(0, 1, 2), dtype=str)
        pol_incident = settings[0].astype(float)
        pol_scattered = settings[1].astype(float)
        axis = settings[2, 0].lower()
    else:
        print("Please define the experimental geometry.")
        pol_incident_str = input("Enter polarization of incident light (e.g., 1.0 0.0 0.0): ")
        pol_scattered_str = input("Enter polarization of scattered light (e.g., 1.0 0.0 0.0): ")
        axis = input("Enter surface normal direction (x, y, or z): ").lower()

        pol_incident = np.array(pol_incident_str.split(), dtype=float)
        pol_scattered = np.array(pol_scattered_str.split(), dtype=float)
        
        with open("input", "w") as f:
            f.write(f"{pol_incident[0]:4.1f} {pol_incident[1]:4.1f} {pol_incident[2]:4.1f}   ! Incident polarization\n")
            f.write(f"{pol_scattered[0]:4.1f} {pol_scattered[1]:4.1f} {pol_scattered[2]:4.1f}   ! Scattered polarization\n")
            f.write(f"{axis} 0.0 0.0      ! Surface normal\n")
            
    return pol_incident, pol_scattered, axis

def read_band_yaml(filepath="band.yaml"):
    """Parses a Phonopy band.yaml file to get phonon modes."""
    print(f"Reading phonon modes from {filepath}...")
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    total_atoms = data['natom']
    n_modes = total_atoms * 3
    
    # Assumes Gamma point calculation is the first entry in the 'phonon' list
    phonons = data['phonon'][0]['band'] 
    
    frequencies = np.array([p['frequency'] for p in phonons])
    eigenvectors_raw = np.array([p['eigenvector'] for p in phonons])
    
    masses = np.array([atom['mass'] for atom in data['points']])

    # Eigenvectors are complex, take the real part and reshape
    eigenvectors = eigenvectors_raw[:, :, :, 0].reshape(n_modes, n_modes)
    
    # Normalize by mass to get eigendisplacements
    masses_expanded = np.repeat(masses, 3) 
    eigendisplacements = eigenvectors / np.sqrt(masses_expanded)
    eigendisplacements = eigendisplacements.reshape(n_modes, total_atoms, 3)
    
    return frequencies, eigendisplacements, masses, total_atoms


def read_dielectric_derivatives(filepath, total_atoms):
    """Parses the dielectric_derivatives to get the atomic Raman tensors (derivatives)."""
    print(f"Reading atomic Raman tensors from {filepath}...")
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    real_start_idx = -1
    imag_start_idx = -1
    for i, line in enumerate(lines):
        if "! The Real Part of Raman tensor:" in line:
            real_start_idx = i + 1
        if "! The Imaginary Part of Raman tensor:" in line:
            imag_start_idx = i + 1

    # Helper to check if a string can be a float
    def is_float(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

def read_dielectric_derivatives(filepath, total_atoms):
    """Parses the dielectric_derivatives to get the atomic Raman tensors (derivatives)."""
    print(f"Reading atomic Raman tensors from {filepath}...")
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    real_start_idx = -1
    imag_start_idx = -1
    for i, line in enumerate(lines):
        if "! The Real Part of dielectric tensor derivatives:" in line:
            real_start_idx = i + 1
        if "! The Imaginary Part of dielectric tensor derivatives:" in line:
            imag_start_idx = i + 1

    def is_float(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def parse_tensor_block(start_idx, end_pattern=None):
        if start_idx == -1: return {}
        
        block_lines = lines[start_idx:]
        tensors = {}
        atom_counter = 0
        
        line_idx = 0
        while line_idx < len(block_lines):
            line = block_lines[line_idx]

            if end_pattern and end_pattern in line:
                break # Stop parsing when we hit the end pattern
            
            parts = line.split()
            
            # Check if line is a potential atom header: not empty, not a comment, starts with non-number
            if parts and "!" not in line and not is_float(parts[0]):
                atom_counter += 1
                atom_idx = atom_counter
                
                # The next 3 lines are the 3x9 tensor data
                tensor_lines = [block_lines[line_idx+1], block_lines[line_idx+2], block_lines[line_idx+3]]
                tensor_3x9 = np.array([list(map(float, l.split())) for l in tensor_lines])
                
                # Correctly slice the 3x9 matrix into three 3x3 tensors for d/dx, d/dy, d/dz
                # and stack them along a new first axis (axis=0)
                tensor_3x3x3 = np.array([
                    tensor_3x9[:, 0:3], # Derivative w.r.t. X displacement
                    tensor_3x9[:, 3:6], # Derivative w.r.t. Y displacement
                    tensor_3x9[:, 6:9]  # Derivative w.r.t. Z displacement
                ])

                tensors[atom_idx] = tensor_3x3x3
                line_idx += 3 # Advance index past the data lines we just read
            
            line_idx += 1
        return tensors

    eps_der_real_dict = parse_tensor_block(real_start_idx, end_pattern="! The Imaginary Part")
    eps_der_imag_dict = parse_tensor_block(imag_start_idx)
    
    symmetry_used = len(eps_der_real_dict) < total_atoms
    
    eps_der_real = np.zeros((total_atoms, 3, 3, 3)) # Shape: (atom, alpha, i, j)
    eps_der_imag = np.zeros((total_atoms, 3, 3, 3))
    
    # Populate the main arrays from the parsed dictionaries (1-based index)
    for i in range(total_atoms):
        if (i + 1) in eps_der_real_dict:
            eps_der_real[i] = eps_der_real_dict.get(i + 1, 0)
            eps_der_imag[i] = eps_der_imag_dict.get(i + 1, 0)

    return eps_der_real, eps_der_imag, symmetry_used

def read_irreps(filepath="irreps.yaml"):
    """Reads irreducible representations from Phonopy's irreps.yaml."""
    if not os.path.exists(filepath):
        return None
    print(f"Reading irreducible representations from {filepath}...")
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    reps = {}
    for mode in data['normal_modes']:
        for band_index in mode['band_indices']:
            reps[band_index] = mode['ir_label']
    return [reps.get(i + 1, '---') for i in range(len(reps))]
    
def run_raman_tensor():
    """Main function to synthesize all data and calculate final intensities."""
    # --- Check for required input files ---
    # Smartly find the dielectric_derivatives based on a pattern
    dielectric_derivatives_path = None
    for f in os.listdir('.'):
        if f.startswith("dielectric_derivatives_"):
            dielectric_derivatives_path = f
            break
    if dielectric_derivatives_path is None:
        print("***** dielectric_derivatives_<freq> not found. Did you run calculate_dielectric_derivatives.py? *****")
        sys.exit(1)

    # --- Get Inputs ---
    pol_incident, pol_scattered, axis = get_user_input()
    frequencies, eigendisps, masses, total_atoms = read_band_yaml()
    eps_der_real, eps_der_imag, symmetry_used = read_dielectric_derivatives(dielectric_derivatives_path, total_atoms)
    representations = read_irreps()
    
    n_modes = total_atoms * 3
    if symmetry_used:
        print("Symmetry was used. Full tensor reconstruction is not yet implemented in this script.")
        print("Results will only be correct if dielectric_derivatives contains data for all atoms.")
        # NOTE: A full translation would require parsing 'symmetry_operation_matrices'
        # and applying the M_inv * T * M transformation, which is complex.
        # This script proceeds assuming the full tensor is available.

    # --- Calculate Mode Raman Tensors ---
    print("Calculating Raman tensors for each mode...")
    # m: mode, j: atom, a: alpha (xyz), i,k: tensor components
    raman_tensor_real = np.einsum('jaik,mja->mik', eps_der_real, eigendisps)
    raman_tensor_imag = np.einsum('jaik,mja->mik', eps_der_imag, eigendisps)
    raman_tensor_cmplx = raman_tensor_real + 1j * raman_tensor_imag

    # --- Calculate Intensities ---
    print("Calculating Raman intensities...")
    # Specific polarization geometry
    # Intensity = | e_s . R_m . e_i |^2
    contracted_tensor = np.einsum('i,mik,k->m', pol_scattered, raman_tensor_cmplx, pol_incident)
    intensities = np.abs(contracted_tensor)**2

    # Polarization-averaged for backscattering
    if axis == 'z':
        avg_intensities = (np.abs(raman_tensor_cmplx[:, 0, 0])**2 + np.abs(raman_tensor_cmplx[:, 0, 1])**2 +
                           np.abs(raman_tensor_cmplx[:, 1, 0])**2 + np.abs(raman_tensor_cmplx[:, 1, 1])**2)
    elif axis == 'y':
        # ... similar logic for other axes ...
        avg_intensities = (np.abs(raman_tensor_cmplx[:, 0, 0])**2 + np.abs(raman_tensor_cmplx[:, 0, 2])**2 +
                           np.abs(raman_tensor_cmplx[:, 2, 0])**2 + np.abs(raman_tensor_cmplx[:, 2, 2])**2)
    elif axis == 'x':
        avg_intensities = (np.abs(raman_tensor_cmplx[:, 1, 1])**2 + np.abs(raman_tensor_cmplx[:, 1, 2])**2 +
                           np.abs(raman_tensor_cmplx[:, 2, 1])**2 + np.abs(raman_tensor_cmplx[:, 2, 2])**2)
    else: # Polycrystalline average
        avg_intensities = np.sum(np.abs(raman_tensor_cmplx)**2, axis=(1, 2))

    # Apply temperature correction (Bose-Einstein factor)
    # h*cm-1/k_B*T at 298K
    const = 0.004824125
    thz_to_cm1 = 33.35641
    freq_cm1 = frequencies * thz_to_cm1
    
    with np.errstate(divide='ignore', invalid='ignore'):
        occupation = 1.0 / (np.exp(freq_cm1 * const) - 1.0)
    occupation[~np.isfinite(occupation)] = 0 # Handle div by zero for low freq
    
    temp_factor = (occupation + 1.0) / frequencies
    temp_factor[frequencies < 0.03] = 1.0 # Avoid division by zero for acoustic modes

    intensities *= temp_factor
    avg_intensities *= temp_factor

    # --- Write Output Files ---
    print("Writing final output files...")
    
    # 1. Raman_tensor file
    with open("Raman_tensor.dat", "w") as f:
        f.write("# Mode   Freq(THz)   Freq(cm-1)   Irrep.   Raman Tensor (Real + i*Imaginary)\n")
        f.write("#--------------------------------------------------------------------------\n")
        for i in range(n_modes):
            rep = representations[i] if representations else "---"
            f.write(f"{i+1:5d} {frequencies[i]:10.3f} {freq_cm1[i]:11.3f}   {rep:<8s}\n")
            for j in range(3):
                row_str = "  ".join([f"{raman_tensor_cmplx[i, j, k].real:10.3f}{raman_tensor_cmplx[i, j, k].imag:+10.3f}j" for k in range(3)])
                f.write(f"    {row_str}\n")
            f.write("\n")

    # 2. Specific intensity file
    with open("Raman_intensity_specific.dat", "w") as f:
        f.write("# Freq(cm-1)   Intensity(arb.)   Irrep.\n")
        for i in range(n_modes):
            rep = representations[i] if representations else ""
            f.write(f"{freq_cm1[i]:12.3f} {intensities[i]:18.4f}   {rep}\n")
    
    # 3. Averaged intensity file
    with open("Raman_intensity_averaged.dat", "w") as f:
        f.write("# Freq(cm-1)   Intensity(arb.)   Irrep.\n")
        for i in range(n_modes):
            rep = representations[i] if representations else ""
            f.write(f"{freq_cm1[i]:12.3f} {avg_intensities[i]:18.4f}   {rep}\n")
            
    print("\ncalculate_spectrum.py finished successfully.")
    print("Generated Raman_tensor.dat, Raman_intensity_specific.dat, and Raman_intensity_averaged.dat")
    print("You can now plot the .dat files to see the spectrum.")

if __name__ == "__main__":
    run_raman_tensor()
