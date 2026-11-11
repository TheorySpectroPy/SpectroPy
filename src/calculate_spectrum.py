import numpy as np
import os
import sys
import re
from spectropy_config import DEFAULT_INPUT, Settings, read_settings, write_settings
from spectropy_phonons import read_registered_modes, read_irreps as _read_irreps


def read_polarization_mode(input_path=DEFAULT_INPUT):
    return read_settings(input_path).polarization


def get_user_input(input_path=DEFAULT_INPUT):
    settings = read_settings(input_path)
    if os.path.exists(input_path):
        print(f"Reading experimental geometry from '{input_path}'...")
        if (
            settings.incident_polarization is None
            or settings.scattered_polarization is None
            or settings.surface_normal is None
        ):
            raise ValueError(f"{input_path} needs incident/scattered polarizations and a surface normal")
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

        write_settings(input_path, Settings(tuple(pol_incident), tuple(pol_scattered), axis))

    return pol_incident, pol_scattered, axis

def read_band_yaml(filepath="band.yaml"):
    print(f"Reading phonon modes from {filepath}...")
    modes = read_registered_modes(filepath, provenance="atom_registration.json")
    return modes.frequencies_thz, modes.eigendisplacements, modes.natoms

def read_dielectric_derivatives(filepath, total_atoms):
    print(f"Reading atomic Raman tensors from {filepath}...")
    if str(filepath).endswith('.json'):
        from spectropy_directional import reconstruct_directional_file
        return reconstruct_directional_file(filepath, total_atoms)
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

            if end_pattern and "Imaginary Part" in line:
                break
            
            parts = line.split()
            
            if parts and "!" not in line and not is_float(parts[0]):
                atom_counter += 1
                label = re.match(r"\s*Atom\s*(\d+)\b", line, re.IGNORECASE)
                atom_idx = int(label.group(1)) if label else atom_counter
                if atom_idx in tensors or not 1 <= atom_idx <= total_atoms:
                    raise ValueError(f"Invalid or duplicate derivative atom index {atom_idx}")
                
                tensor_lines = [block_lines[line_idx+1], block_lines[line_idx+2], block_lines[line_idx+3]]
                tensor_3x9 = np.array([list(map(float, l.split())) for l in tensor_lines])
                if tensor_3x9.shape != (3, 9) or not np.isfinite(tensor_3x9).all():
                    raise ValueError(f'Malformed or nonfinite derivative tensor for atom {atom_idx}')
                
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
    
    missing = [i + 1 for i in range(total_atoms)
               if (i + 1) not in eps_der_real_dict or (i + 1) not in eps_der_imag_dict]
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
    
def run_raman_tensor(dielectric_derivatives_path=None, _phonons=None, output_prefix="", input_path=DEFAULT_INPUT):
    if dielectric_derivatives_path is None:
        for f in sorted(os.listdir('.')):
            if f.startswith("epsilon_derivative_"):
                dielectric_derivatives_path = f
                break
    if dielectric_derivatives_path is None:
        print("***** epsilon_derivative_<freq> not found. Did you run `spectropy derivatives`? *****")
        sys.exit(1)

    pol_incident, pol_scattered, axis = get_user_input(input_path)
    polarization_mode = read_polarization_mode(input_path)
    if _phonons is None:
        frequencies, eigendisps, total_atoms = read_band_yaml()
    else:
        frequencies, eigendisps, total_atoms = (_phonons.frequencies_thz,
            _phonons.eigendisplacements, _phonons.natoms)
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
    cutoff_cm1 = read_settings(input_path).frequency_cutoff_cm1
    if cutoff_cm1 is not None:
        selected = (frequencies > 0.0) & (frequencies >= cutoff_cm1 / thz_to_cm1)
        occupation = 1.0 / (np.exp(freq_cm1[selected] * 0.004824125) - 1.0)
        factor = (occupation + 1.0) / frequencies[selected]
        intensities[selected] *= factor
        if avg_intensities is not None:
            avg_intensities[selected] *= factor

    print("Writing final output files...")
    
    filename = os.path.basename(dielectric_derivatives_path)
    energy = (filename.removeprefix('raw_directional_derivative_').removesuffix('.json')
              if filename.startswith('raw_directional_derivative_')
              else filename.split('epsilon_derivative_', 1)[1])

    def representation(index):
        if representations is not None and index < len(representations):
            return representations[index]
        return "---"

    tensor_path = f"{output_prefix}Raman_tensor_{energy}eV"
    with open(tensor_path, "w") as f:
        f.write("# Mode   Freq(THz)   Freq(cm-1)   Irrep.   Raman Tensor (Real + i*Imaginary)\n")
        f.write("#--------------------------------------------------------------------------\n")
        for i in range(n_modes):
            rep = representation(i)
            f.write(f"{i+1:5d} {frequencies[i]:10.3f} {freq_cm1[i]:11.3f}   {rep:<8s}\n")
            for j in range(3):
                row_str = "  ".join([f"{raman_tensor_cmplx[i, j, k].real:10.3f}{raman_tensor_cmplx[i, j, k].imag:+10.3f}j" for k in range(3)])
                f.write(f"    {row_str}\n")
            f.write("\n")

    # Preserve the legacy active tensor alias, without overwriting named variants.
    if not output_prefix:
        import shutil
        shutil.copyfile(tensor_path, "Raman_tensor")

    with open(f"{output_prefix}Raman_intensity_complex_{energy}eV", "w") as f:
        for i in range(n_modes):
            f.write(f"{freq_cm1[i]:12.3f} {intensities[i]:18.4f}       {representation(i)}\n")

    if avg_intensities is not None:
        with open(f"{output_prefix}Raman_intensity_polarization_averaged_{energy}eV", "w") as f:
            for i in range(n_modes):
                f.write(f"{freq_cm1[i]:12.3f} {avg_intensities[i]:18.4f}       {representation(i)}\n")
            
    print("\ncalculate_spectrum.py finished successfully.")
    print(f"Generated Raman_tensor and intensity files for {energy} eV.")


def run_raman_tensors_for_input(input_path=DEFAULT_INPUT, _phonons=None):
    from spectropy_config import read_settings

    settings = read_settings(input_path)
    for energy in settings.laser_energies:
        derivative_path = f"epsilon_derivative_{energy:.2f}"
        raw_path = f"raw_{derivative_path}"
        directional_path = f"raw_directional_derivative_{energy:.2f}.json"
        if os.path.isfile(directional_path):
            raw_path = directional_path
        corrected_path = f"corrected_{derivative_path}"
        # Older unprojected datasets may only have the legacy filename. Never
        # label an active, potentially corrected field as raw when denoising.
        if not os.path.isfile(raw_path) and not settings.site_symmetry_projection:
            raw_path = derivative_path
        required = [raw_path]
        if settings.site_symmetry_projection:
            required.append(corrected_path)
        for path in required:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Missing {path}; run `spectropy derivatives` first")
        run_raman_tensor(raw_path, _phonons=_phonons, output_prefix="raw_", input_path=input_path)
        if settings.site_symmetry_projection:
            run_raman_tensor(corrected_path, _phonons=_phonons, output_prefix="corrected_", input_path=input_path)
        # Existing workflow consumers continue to see the selected field.
        run_raman_tensor(corrected_path if settings.site_symmetry_projection else raw_path,
                         _phonons=_phonons, input_path=input_path)

if __name__ == "__main__":
    run_raman_tensors_for_input()
