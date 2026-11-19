import numpy as np
import os
import sys
import re
from .config import DEFAULT_INPUT, resolve_settings
from .phonons import read_registered_modes, read_irreps as _read_irreps


def polarization_vectors(settings):
    """Validate configured geometry without prompting or modifying input.yaml."""
    if settings.polarization != 'specific':
        raise ValueError('polarization.mode must be specific')
    if settings.incident_polarization is None or settings.scattered_polarization is None:
        raise ValueError('polarization.incident and polarization.scattered must be supplied in input.yaml')
    return (np.asarray(settings.incident_polarization),
            np.asarray(settings.scattered_polarization), settings.surface_normal)


def read_band_yaml(filepath="band.yaml"):
    print(f"Reading phonon modes from {filepath}...")
    modes = read_registered_modes(filepath)
    return modes.frequencies_thz, modes.eigendisplacements, modes.natoms

def read_dielectric_derivatives(filepath, total_atoms):
    print(f"Reading Cartesian atomic dielectric derivatives from {filepath}...")
    with open(filepath, 'r') as f:
        lines = f.readlines()
    if lines and lines[0].startswith('! Raw directional dielectric derivatives'):
        raise ValueError('Raw measured-direction derivatives are not a Cartesian atomic field; '
                         'use the final epsilon_derivative_<E> file for Cartesian contraction')
        
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
              "regenerate the final atomic derivatives with `spectropy tensors`. *****")
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
    
def run_raman_tensor(dielectric_derivatives_path, _phonons=None, output_prefix="", input_path=DEFAULT_INPUT):
    """Import an existing atomic derivative file (Python helper for bundled data).

    Files already contain Vcell/(4 pi); do not apply this factor again.
    The normal CLI workflow contracts the reconstructed field in memory.
    """
    from .tensors import write_raman_tensors

    if _phonons is None:
        frequencies, eigendisps, total_atoms = read_band_yaml()
    else:
        frequencies, eigendisps, total_atoms = (_phonons.frequencies_thz,
            _phonons.eigendisplacements, _phonons.natoms)
    real, imaginary = read_dielectric_derivatives(dielectric_derivatives_path, total_atoms)
    tensors = np.einsum('jaik,mja->mik', real + 1j * imaginary, eigendisps)
    name = os.path.basename(dielectric_derivatives_path)
    energy = float(name.rsplit('_', 1)[-1].removesuffix('.json'))
    write_raman_tensors(f"{output_prefix}Raman_tensor_{energy:.2f}eV",
                        frequencies, tensors, read_irreps())


def calculate_intensities(frequencies, tensors, settings):
    """Return frequencies in cm^-1 and intensities, using one shared convention."""
    from .tensors import THZ_TO_CM1

    pol_incident, pol_scattered, _ = polarization_vectors(settings)
    intensities = np.abs(np.einsum('i,mik,k->m', pol_scattered, tensors, pol_incident))**2
    freq_cm1 = frequencies * THZ_TO_CM1
    cutoff = settings.frequency_cutoff_cm1
    factor = np.ones_like(frequencies, dtype=float)
    if cutoff is not None:
        selected = (frequencies > 0.0) & (freq_cm1 >= cutoff)
        factor[selected] /= frequencies[selected]
    if settings.temperature > 0:
        positive = frequencies > 0.0
        # hc/k_B in cm K; stable Stokes n+1, including the low-T limit.
        exponent = freq_cm1[positive] * 1.4387768775039338 / settings.temperature
        factor[positive] /= -np.expm1(-exponent)
    intensities *= factor
    return freq_cm1, {"complex": intensities}


def write_intensities(energy, frequencies, tensors, labels, settings):
    freq_cm1, outputs = calculate_intensities(frequencies, tensors, settings)
    for kind, values in outputs.items():
        path = f"Raman_intensity_{kind}_{energy:.2f}eV"
        with open(path, "w") as stream:
            for frequency, intensity, label in zip(freq_cm1, values, labels):
                stream.write(f"{frequency:.16e} {intensity:.16e} {label}\n")
        print(f"Wrote {path}")


def run_spectrum(input_path=DEFAULT_INPUT, *, settings=None):
    """Calculate intensities from final mode tensors, not atomic intermediates."""
    from .tensors import read_raman_tensors

    settings = resolve_settings(settings, input_path)
    paths = [f"Raman_tensor_{energy:.2f}eV" for energy in settings.laser_energies]
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {path}; run `spectropy tensors` first")
    for energy, path in zip(settings.laser_energies, paths):
        mode_indices, frequencies, tensors, labels = read_raman_tensors(path)
        write_intensities(energy, frequencies, tensors, labels, settings)


if __name__ == "__main__":
    run_spectrum()
