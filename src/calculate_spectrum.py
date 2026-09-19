import numpy as np
import os
import sys
from spectropy_config import Settings, read_settings, write_settings
from spectropy_phonons import read_gamma_modes, read_irreps as _read_irreps


def read_polarization_mode(input_path="input"):
    return read_settings(input_path).polarization


def get_user_input():
    if os.path.exists("input"):
        print("Reading experimental geometry from 'input' file...")
        settings = read_settings("input")
        if (
            settings.incident_polarization is None
            or settings.scattered_polarization is None
            or settings.surface_normal is None
        ):
            raise ValueError("input needs incident/scattered polarizations and a surface normal")
        pol_incident = np.array(settings.incident_polarization)
        pol_scattered = np.array(settings.scattered_polarization)
        axis = settings.surface_normal
    else:
        print("Please define the experimental geometry.")
        pol_incident_str = input("Enter polarization of incident light (e.g., 1.0 0.0 0.0): ")
        pol_scattered_str = input("Enter polarization of scattered light (e.g., 1.0 0.0 0.0): ")
        axis = input("Enter surface normal direction (x, y, or z): ").lower()

        pol_incident = np.array(pol_incident_str.split(), dtype=float)
        pol_scattered = np.array(pol_scattered_str.split(), dtype=float)

        write_settings("input", Settings(tuple(pol_incident), tuple(pol_scattered), axis))

    return pol_incident, pol_scattered, axis

def read_band_yaml(filepath="band.yaml"):
    print(f"Reading phonon modes from {filepath}...")
    modes = read_gamma_modes(filepath)
    return modes.frequencies_thz, modes.eigendisplacements, modes.natoms

def read_dielectric_derivatives(filepath, total_atoms):
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
                break
            
            parts = line.split()
            
            if parts and "!" not in line and not is_float(parts[0]):
                atom_counter += 1
                atom_idx = atom_counter
                
                tensor_lines = [block_lines[line_idx+1], block_lines[line_idx+2], block_lines[line_idx+3]]
                tensor_3x9 = np.array([list(map(float, l.split())) for l in tensor_lines])
                
                tensor_3x3x3 = np.array([
                    tensor_3x9[:, 0:3],
                    tensor_3x9[:, 3:6],
                    tensor_3x9[:, 6:9]
                ])

                tensors[atom_idx] = tensor_3x3x3
                line_idx += 3
            line_idx += 1
        return tensors

    eps_der_real_dict = parse_tensor_block(real_start_idx, end_pattern="! The Imaginary Part")
    eps_der_imag_dict = parse_tensor_block(imag_start_idx)
    
    missing = [i + 1 for i in range(total_atoms) if (i + 1) not in eps_der_real_dict]
    if missing:
        print(f"***** dielectric_derivatives file is missing atom(s) {missing} -- "
              "was it written by an up-to-date calculate_dielectric_derivatives.py? *****")
        sys.exit(1)

    eps_der_real = np.zeros((total_atoms, 3, 3, 3))
    eps_der_imag = np.zeros((total_atoms, 3, 3, 3))

    for i in range(total_atoms):
        eps_der_real[i] = eps_der_real_dict[i + 1]
        eps_der_imag[i] = eps_der_imag_dict[i + 1]

    return eps_der_real, eps_der_imag

def read_irreps(filepath="irreps.yaml"):
    labels = _read_irreps(filepath)
    if labels is not None:
        print(f"Reading irreducible representations from {filepath}...")
    return labels
    
def run_raman_tensor(dielectric_derivatives_path=None):
    if dielectric_derivatives_path is None:
        for f in sorted(os.listdir('.')):
            if f.startswith("epsilon_derivative_"):
                dielectric_derivatives_path = f
                break
    if dielectric_derivatives_path is None:
        print("***** epsilon_derivative_<freq> not found. Did you run `spectropy derivatives`? *****")
        sys.exit(1)

    pol_incident, pol_scattered, axis = get_user_input()
    polarization_mode = read_polarization_mode()
    frequencies, eigendisps, total_atoms = read_band_yaml()
    eps_der_real, eps_der_imag = read_dielectric_derivatives(dielectric_derivatives_path, total_atoms)
    representations = read_irreps()

    n_modes = total_atoms * 3

    print("Calculating Raman tensors for each mode...")
    raman_tensor_real = np.einsum('jaik,mja->mik', eps_der_real, eigendisps)
    raman_tensor_imag = np.einsum('jaik,mja->mik', eps_der_imag, eigendisps)
    raman_tensor_cmplx = raman_tensor_real + 1j * raman_tensor_imag

    print("Calculating Raman intensities...")
    contracted_tensor = np.einsum('i,mik,k->m', pol_scattered, raman_tensor_cmplx, pol_incident)
    intensities = np.abs(contracted_tensor)**2

    avg_intensities = None
    if polarization_mode == "average":
        if axis == 'z':
            avg_intensities = (np.abs(raman_tensor_cmplx[:, 0, 0])**2 + np.abs(raman_tensor_cmplx[:, 0, 1])**2 +
                               np.abs(raman_tensor_cmplx[:, 1, 0])**2 + np.abs(raman_tensor_cmplx[:, 1, 1])**2)
        elif axis == 'y':
            avg_intensities = (np.abs(raman_tensor_cmplx[:, 0, 0])**2 + np.abs(raman_tensor_cmplx[:, 0, 2])**2 +
                               np.abs(raman_tensor_cmplx[:, 2, 0])**2 + np.abs(raman_tensor_cmplx[:, 2, 2])**2)
        elif axis == 'x':
            avg_intensities = (np.abs(raman_tensor_cmplx[:, 1, 1])**2 + np.abs(raman_tensor_cmplx[:, 1, 2])**2 +
                               np.abs(raman_tensor_cmplx[:, 2, 1])**2 + np.abs(raman_tensor_cmplx[:, 2, 2])**2)
        else:
            avg_intensities = np.sum(np.abs(raman_tensor_cmplx)**2, axis=(1, 2))

    thz_to_cm1 = 33.35641
    freq_cm1 = frequencies * thz_to_cm1
    cutoff_cm1 = read_settings().frequency_cutoff_cm1
    if cutoff_cm1 is not None:
        selected = (frequencies > 0.0) & (frequencies >= cutoff_cm1 / thz_to_cm1)
        occupation = 1.0 / (np.exp(freq_cm1[selected] * 0.004824125) - 1.0)
        factor = (occupation + 1.0) / frequencies[selected]
        intensities[selected] *= factor
        if avg_intensities is not None:
            avg_intensities[selected] *= factor

    print("Writing final output files...")
    
    energy = dielectric_derivatives_path.removeprefix("epsilon_derivative_")

    def representation(index):
        if representations is not None and index < len(representations):
            return representations[index]
        return "---"

    with open("Raman_tensor", "w") as f:
        f.write("# Mode   Freq(THz)   Freq(cm-1)   Irrep.   Raman Tensor (Real + i*Imaginary)\n")
        f.write("#--------------------------------------------------------------------------\n")
        for i in range(n_modes):
            rep = representation(i)
            f.write(f"{i+1:5d} {frequencies[i]:10.3f} {freq_cm1[i]:11.3f}   {rep:<8s}\n")
            for j in range(3):
                row_str = "  ".join([f"{raman_tensor_cmplx[i, j, k].real:10.3f}{raman_tensor_cmplx[i, j, k].imag:+10.3f}j" for k in range(3)])
                f.write(f"    {row_str}\n")
            f.write("\n")

    with open(f"Raman_intensity_complex_{energy}eV", "w") as f:
        for i in range(n_modes):
            f.write(f"{freq_cm1[i]:12.3f} {intensities[i]:18.4f}       {representation(i)}\n")

    if avg_intensities is not None:
        with open(f"Raman_intensity_polarization_averaged_{energy}eV", "w") as f:
            for i in range(n_modes):
                f.write(f"{freq_cm1[i]:12.3f} {avg_intensities[i]:18.4f}       {representation(i)}\n")
            
    print("\ncalculate_spectrum.py finished successfully.")
    print(f"Generated Raman_tensor and intensity files for {energy} eV.")


def run_raman_tensors_for_input(input_path="input"):
    from spectropy_config import read_settings

    for energy in read_settings(input_path).laser_energies:
        derivative_path = f"epsilon_derivative_{energy:.2f}"
        if not os.path.isfile(derivative_path):
            raise FileNotFoundError(f"Missing {derivative_path}; run `spectropy derivatives` first")
        run_raman_tensor(derivative_path)

if __name__ == "__main__":
    run_raman_tensors_for_input()
