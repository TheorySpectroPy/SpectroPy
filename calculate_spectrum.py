import numpy as np
import os
import sys
import yaml
import re

_SETTINGS_KEY_RE = re.compile(
    r"^(laser_energies?|broadening_fwhm|broadening_type|displacement_amplitude|polarization)\s*(?::|=)",
    re.I,
)


def read_polarization_mode(input_path="input"):
    """
    Read ``polarization`` from ``input``: ``average`` requests the
    orientation-averaged (powder) intensity in addition to the specific
    incident/scattered geometry; anything else (the default, "specific")
    means only the specific-polarization intensity is computed -- averaging
    is never done unless explicitly asked for.
    """
    if not os.path.exists(input_path):
        return "specific"
    with open(input_path) as handle:
        for raw_line in handle:
            line = raw_line.split("!", 1)[0].split("#", 1)[0].strip()
            match = re.match(r"^polarization\s*(?::|=)\s*(.+)$", line, re.I)
            if match:
                value = match.group(1).strip().lower()
                if value not in ("average", "specific"):
                    raise ValueError(
                        f"polarization must be 'average' or 'specific', got: {raw_line.rstrip()}"
                    )
                return value
    return "specific"


def get_user_input():
    """Gets polarization vectors from 'input' file or user prompt."""
    if os.path.exists("input"):
        print("Reading experimental geometry from 'input' file...")
        settings = []
        with open("input") as handle:
            for raw_line in handle:
                line = raw_line.split("!", 1)[0].split("#", 1)[0].strip()
                if line and not _SETTINGS_KEY_RE.match(line):
                    settings.append(line.split())
        if len(settings) < 3:
            raise ValueError("input needs incident/scattered polarizations and a surface normal")
        pol_incident = np.array(settings[0][:3], dtype=float)
        pol_scattered = np.array(settings[1][:3], dtype=float)
        axis = settings[2][0].lower()
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
            f.write("laser_energies: 0.00 ! eV\n")
            f.write("broadening_fwhm: 1.0 ! cm-1\n")
            f.write("broadening_type: lorentzian\n")
            f.write("polarization: specific ! 'average' to also compute the orientation-averaged intensity\n")

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
        if "Real Part" in line:
            real_start_idx = i + 1
        if "Imaginary Part" in line:
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
    
    # calculate_dielectric_derivatives.py now always writes data for every
    # atom (symmetry-equivalent atoms not directly displaced are filled in
    # via reconstruct_dielectric_derivatives.expand_to_equivalent_atoms
    # before the file is written), so no separate "symmetry was used, full
    # reconstruction not implemented" fallback path is needed here anymore.
    missing = [i + 1 for i in range(total_atoms) if (i + 1) not in eps_der_real_dict]
    if missing:
        print(f"***** dielectric_derivatives file is missing atom(s) {missing} -- "
              "was it written by an up-to-date calculate_dielectric_derivatives.py? *****")
        sys.exit(1)

    eps_der_real = np.zeros((total_atoms, 3, 3, 3)) # Shape: (atom, alpha, i, j)
    eps_der_imag = np.zeros((total_atoms, 3, 3, 3))

    for i in range(total_atoms):
        eps_der_real[i] = eps_der_real_dict[i + 1]
        eps_der_imag[i] = eps_der_imag_dict[i + 1]

    return eps_der_real, eps_der_imag

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
    
def run_raman_tensor(dielectric_derivatives_path=None):
    """Main function to synthesize all data and calculate final intensities."""
    # --- Check for required input files ---
    # Smartly find the dielectric_derivatives based on a pattern
    if dielectric_derivatives_path is None:
        for f in sorted(os.listdir('.')):
            if f.startswith("epsilon_derivative_"):
                dielectric_derivatives_path = f
                break
    if dielectric_derivatives_path is None:
        print("***** epsilon_derivative_<freq> not found. Did you run `spectropy derivatives`? *****")
        sys.exit(1)

    # --- Get Inputs ---
    pol_incident, pol_scattered, axis = get_user_input()
    polarization_mode = read_polarization_mode()
    frequencies, eigendisps, masses, total_atoms = read_band_yaml()
    eps_der_real, eps_der_imag = read_dielectric_derivatives(dielectric_derivatives_path, total_atoms)
    representations = read_irreps()

    n_modes = total_atoms * 3

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

    # Polarization-averaged for backscattering -- only computed when
    # explicitly requested via "polarization: average" in input (see
    # read_polarization_mode); otherwise only the specific incident/scattered
    # geometry above is used, and no averaged file is written at all.
    avg_intensities = None
    if polarization_mode == "average":
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
    if avg_intensities is not None:
        avg_intensities *= temp_factor

    # --- Write Output Files ---
    print("Writing final output files...")
    
    energy = dielectric_derivatives_path.removeprefix("epsilon_derivative_")

    # Match the established Raman-workflow filenames: a single Raman_tensor
    # file plus one two-column intensity file per laser energy.
    with open("Raman_tensor", "w") as f:
        f.write("# Mode   Freq(THz)   Freq(cm-1)   Irrep.   Raman Tensor (Real + i*Imaginary)\n")
        f.write("#--------------------------------------------------------------------------\n")
        for i in range(n_modes):
            rep = representations[i] if representations else "---"
            f.write(f"{i+1:5d} {frequencies[i]:10.3f} {freq_cm1[i]:11.3f}   {rep:<8s}\n")
            for j in range(3):
                row_str = "  ".join([f"{raman_tensor_cmplx[i, j, k].real:10.3f}{raman_tensor_cmplx[i, j, k].imag:+10.3f}j" for k in range(3)])
                f.write(f"    {row_str}\n")
            f.write("\n")

    with open(f"Raman_intensity_complex_{energy}eV", "w") as f:
        for i in range(n_modes):
            f.write(f"{freq_cm1[i]:12.3f} {intensities[i]:18.4f}\n")

    if avg_intensities is not None:
        with open(f"Raman_intensity_polarization_averaged_{energy}eV", "w") as f:
            for i in range(n_modes):
                f.write(f"{freq_cm1[i]:12.3f} {avg_intensities[i]:18.4f}\n")
            
    print("\ncalculate_spectrum.py finished successfully.")
    print(f"Generated Raman_tensor and intensity files for {energy} eV.")


def run_raman_tensors_for_input(input_path="input"):
    """Generate Raman results for every laser energy requested in ``input``."""
    from calculate_dielectric_derivatives import read_laser_energies

    for energy in read_laser_energies(input_path):
        derivative_path = f"epsilon_derivative_{energy:.2f}"
        if not os.path.isfile(derivative_path):
            raise FileNotFoundError(f"Missing {derivative_path}; run `spectropy derivatives` first")
        run_raman_tensor(derivative_path)

if __name__ == "__main__":
    run_raman_tensors_for_input()
