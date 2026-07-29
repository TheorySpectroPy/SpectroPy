import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import sys
import re
import shutil
from spectropy_config import read_settings

# --- 1. Publication-Style Configuration ---
# rc('text', usetex=True) never raises by itself -- it just sets a param;
# matplotlib only discovers latex is missing later, while actually rendering
# text (deep inside tight_layout/savefig), so a try/except around rc() alone
# never catches the missing-latex case. Check for the binary up front instead.
if shutil.which('latex'):
    rc('text', usetex=True)
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica', 'Arial']})
else:
    print("Notice: LaTeX not found. Using standard Matplotlib fonts.")
    rc('text', usetex=False)
    rc('font', family='sans-serif')

# --- 2. Math & Helper Functions ---

def gaussian(x, center, amplitude, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return amplitude * np.exp(-((x - center)**2) / (2 * sigma**2))

def lorentzian(x, center, amplitude, fwhm):
    hwhm = fwhm / 2.0
    return amplitude * (hwhm**2) / ((x - center)**2 + hwhm**2)


def broaden_spectrum(freqs, intensities, fwhm, b_type):
    """Return an unnormalized broadened spectrum on the plotting grid."""
    x_dense = np.linspace(max(0, min(freqs) - 50), max(freqs) + 50, 2000)
    y_dense = np.zeros_like(x_dense)
    profile = gaussian if b_type == 'g' else lorentzian
    for frequency, intensity in zip(freqs, intensities):
        y_dense += profile(x_dense, frequency, intensity, fwhm)
    return x_dense, y_dense


def write_broadened_spectrum(input_file, fwhm, b_type):
    """Write the Raman-workflow-compatible numerical broadening output."""
    raw_data = np.loadtxt(input_file, dtype=float)
    if raw_data.size == 0:
        return None
    if raw_data.ndim == 1:
        raw_data = raw_data.reshape(1, -1)
    x_dense, y_dense = broaden_spectrum(raw_data[:, 0], raw_data[:, 1], fwhm, b_type)

    basename = os.path.basename(input_file)
    energy = basename.removeprefix("Raman_intensity_complex_")
    output_path = os.path.join(
        os.path.dirname(input_file), f"Raman_intensity_complex_broadening_{energy}"
    )
    np.savetxt(output_path, np.column_stack((x_dense, y_dense)), fmt="%12.6f %18.8e")
    print(f"   -> Created {output_path}")
    return output_path

def format_mode_for_latex(mode_str):
    r"""
    Robust formatter. Handles subscripts and primes better.
    Ex: E2g -> E_{2g},  A' -> A^{\prime}
    """
    # Regex to capture: 1) Base letters, 2) Subscript numbers/letters, 3) Primes
    match = re.match(r"^([A-Za-z]+)((?:\d+[a-zA-Z]*)?)((?:\'|\")*)$", mode_str)

    if match:
        base, subscript, primes = match.groups()
        subscript_latex = f"_{{{subscript}}}" if subscript else ""

        if primes == "'":
            prime_latex = r"^{\prime}"
        elif primes == "''":
            prime_latex = r"^{\prime\prime}"
        else:
            prime_latex = ""

        return f'${base}{subscript_latex}{prime_latex}$'

    # Fallback
    return f'${mode_str}$'


def read_broadening_settings(input_path="input"):
    settings = read_settings(input_path)
    return settings.broadening_fwhm, "g" if settings.broadening_type == "gaussian" else "l"

# --- 3. The Plotting Core ---

def process_and_plot(input_file, fwhm=5.0, b_type='l'):
    try:
        # raman_tensor's Raman_intensity_polarization_averaged_<eV>eV /
        # Raman_intensity_complex_<eV>eV have no header row and only two
        # columns: frequency (cm-1) and intensity. There is no per-mode
        # symmetry label in this output (raman_tensor doesn't write an
        # irreps.yaml), so peaks are annotated with their frequency instead.
        raw_data = np.loadtxt(input_file, dtype=float)
        if raw_data.size == 0: return
        if raw_data.ndim == 1: raw_data = raw_data.reshape(1, -1)

        freqs = raw_data[:, 0]
        intensities = raw_data[:, 1]
        modes = [f"{f:.1f}" for f in freqs]
    except Exception as e:
        print(f"      Error reading {os.path.basename(input_file)}: {e}")
        return

    # Generate the same unnormalized curve that is written for the raw
    # polarization-specific spectrum below, then normalize only the plot.
    x_dense, y_dense = broaden_spectrum(freqs, intensities, fwhm, b_type)

    # Normalize
    if np.max(y_dense) > 0:
        y_dense /= np.max(y_dense)
        intensities /= np.max(intensities)

    # --- Plot Setup ---
    # Standard single-column width
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    
    # Professional Blue color
    line_color = '#2c7bb6' 
    
    # Plot Curve
    ax.plot(x_dense, y_dense, color=line_color, linewidth=1.5)
    ax.fill_between(x_dense, y_dense, color=line_color, alpha=0.1)

    # --- Annotations with Arrows ---
    for f, i, m in zip(freqs, intensities, modes):
        # Only label peaks > 10% intensity
        if i > 0.1: 
            y_curve = np.interp(f, x_dense, y_dense)
            formatted_mode = format_mode_for_latex(m)
            
            # Annotation Style
            ax.annotate(formatted_mode,
                        xy=(f, y_curve), 
                        xytext=(f, y_curve + 0.15),
                        fontsize=10,
                        ha='center',
                        arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=3, headlength=3))

    # --- Axis Styling ---
    ax.set_xlabel(r'Raman Shift (cm$^{-1}$)', fontsize=11)
    ax.set_ylabel(r'Intensity (Arb. Units)', fontsize=11)
    
    # Inward ticks
    ax.tick_params(axis='both', which='major', labelsize=10, direction='in', top=True, right=True)
    
    # Hide Y-axis numbers
    ax.set_yticks([])
    
    # Vertical Dashed Grid
    ax.grid(visible=True, which='major', axis='x', linestyle='--', linewidth=0.5, alpha=0.7)
    
    # Set limits
    ax.set_ylim(bottom=-0.02, top=1.35)
    ax.set_xlim(left=0, right=max(freqs)+60)

    plt.tight_layout()
    
    # Save Output -- one file per input (a directory can hold one plot per
    # incident energy), named after the source data file. Not splitext: names
    # like "..._2.33eV" have a '.' inside the energy value itself, which
    # splitext would mistake for a file extension and truncate.
    base = os.path.basename(input_file)
    out_name = os.path.join(os.path.dirname(input_file), f"Raman_plot_styled_{base}.png")
    plt.savefig(out_name, dpi=300)
    plt.close()
    print(f"   -> Created {out_name}")

# --- 4. Automation Logic ---

def run_automation():
    base_path = os.getcwd()
    print(f"--- Automated Raman Plotter (Style: Publication) ---")
    print(f"Scanning: {base_path}")
    
    fwhm, b_type = read_broadening_settings()
    print(f"Broadening: {'Lorentzian' if b_type == 'l' else 'Gaussian'}, FWHM = {fwhm:g} cm-1")

    count = 0
    broadened_count = 0
    # Walk Directory -- raman_tensor names its per-energy output
    # Raman_intensity_complex_<eV>eV (the specific incident/scattered
    # geometry from input, always written) and, only when "polarization:
    # average" was requested (see calculate_spectrum.read_polarization_mode),
    # Raman_intensity_polarization_averaged_<eV>eV too. Plot whichever of the
    # two actually exists for a given energy -- average is only present when
    # explicitly asked for, so this never silently substitutes one for the
    # other.
    avg_re = re.compile(r"^Raman_intensity_polarization_averaged_.+eV$")
    raw_complex_re = re.compile(r"^Raman_intensity_complex_(?!broadening_).+eV$")
    for root, dirs, files in os.walk(base_path):
        for fname in files:
            if raw_complex_re.match(fname):
                full_path = os.path.join(root, fname)
                try:
                    write_broadened_spectrum(full_path, fwhm, b_type)
                    broadened_count += 1
                except Exception as error:
                    print(f"      Error broadening {fname}: {error}")
                print(f"Processing: {os.path.join(os.path.basename(root), fname)}")
                process_and_plot(full_path, fwhm, b_type)
                count += 1
            if avg_re.match(fname):
                full_path = os.path.join(root, fname)
                print(f"Processing: {os.path.join(os.path.basename(root), fname)}")
                process_and_plot(full_path, fwhm, b_type)
                count += 1

    print(f"\nSuccess! Generated {count} plots and {broadened_count} numerical broadened spectra.")

if __name__ == "__main__":
    run_automation()
