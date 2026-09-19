import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import re
import shutil
from spectropy_config import read_settings

if shutil.which('latex'):
    rc('text', usetex=True)
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica', 'Arial']})
else:
    print("Notice: LaTeX not found. Using standard Matplotlib fonts.")
    rc('text', usetex=False)
    rc('font', family='sans-serif')

def gaussian(x, center, amplitude, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return amplitude * np.exp(-((x - center)**2) / (2 * sigma**2))

def lorentzian(x, center, amplitude, fwhm):
    hwhm = fwhm / 2.0
    return amplitude * (hwhm**2) / ((x - center)**2 + hwhm**2)


def broaden_spectrum(freqs, intensities, fwhm, b_type):
    x_dense = np.linspace(max(0, min(freqs) - 50), max(freqs) + 50, 2000)
    y_dense = np.zeros_like(x_dense)
    profile = gaussian if b_type == 'g' else lorentzian
    for frequency, intensity in zip(freqs, intensities):
        y_dense += profile(x_dense, frequency, intensity, fwhm)
    return x_dense, y_dense


def read_intensity_file(input_file):
    frequencies = []
    intensities = []
    modes = []
    with open(input_file) as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) < 2:
                raise ValueError(f"line {line_number} has fewer than two columns")
            frequencies.append(float(fields[0]))
            intensities.append(float(fields[1]))
            modes.append(fields[2] if len(fields) >= 3 else f"{float(fields[0]):.1f}")
    return np.asarray(frequencies), np.asarray(intensities), modes


def write_broadened_spectrum(input_file, fwhm, b_type):
    frequencies, intensities, _ = read_intensity_file(input_file)
    if frequencies.size == 0:
        return None
    x_dense, y_dense = broaden_spectrum(frequencies, intensities, fwhm, b_type)

    basename = os.path.basename(input_file)
    energy = basename.removeprefix("Raman_intensity_complex_")
    output_path = os.path.join(
        os.path.dirname(input_file), f"Raman_intensity_complex_broadening_{energy}"
    )
    np.savetxt(output_path, np.column_stack((x_dense, y_dense)), fmt="%12.6f %18.8e")
    print(f"   -> Created {output_path}")
    return output_path

def format_mode_for_latex(mode_str):
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

    return f'${mode_str}$'


def read_broadening_settings(input_path="input"):
    settings = read_settings(input_path)
    return settings.broadening_fwhm, "g" if settings.broadening_type == "gaussian" else "l"

def process_and_plot(input_file, fwhm=5.0, b_type='l'):
    try:
        freqs, intensities, modes = read_intensity_file(input_file)
        if freqs.size == 0:
            return False
    except Exception as e:
        print(f"      Error reading {os.path.basename(input_file)}: {e}")
        return False

    x_dense, y_dense = broaden_spectrum(freqs, intensities, fwhm, b_type)

    if np.max(y_dense) > 0:
        y_dense /= np.max(y_dense)
        intensities /= np.max(intensities)

    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    
    line_color = '#2c7bb6' 
    
    ax.plot(x_dense, y_dense, color=line_color, linewidth=1.5)
    ax.fill_between(x_dense, y_dense, color=line_color, alpha=0.1)

    for f, i, m in zip(freqs, intensities, modes):
        if i > 0.1: 
            y_curve = np.interp(f, x_dense, y_dense)
            formatted_mode = format_mode_for_latex(m)
            
            ax.annotate(formatted_mode,
                        xy=(f, y_curve), 
                        xytext=(f, y_curve + 0.15),
                        fontsize=10,
                        ha='center',
                        arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=3, headlength=3))

    ax.set_xlabel(r'Raman Shift (cm$^{-1}$)', fontsize=11)
    ax.set_ylabel(r'Intensity (Arb. Units)', fontsize=11)
    
    ax.tick_params(axis='both', which='major', labelsize=10, direction='in', top=True, right=True)
    
    ax.set_yticks([])
    
    ax.grid(visible=True, which='major', axis='x', linestyle='--', linewidth=0.5, alpha=0.7)
    
    ax.set_ylim(bottom=-0.02, top=1.35)
    ax.set_xlim(left=0, right=max(freqs)+60)

    plt.tight_layout()
    
    base = os.path.basename(input_file)
    out_name = os.path.join(os.path.dirname(input_file), f"Raman_plot_styled_{base}.png")
    plt.savefig(out_name, dpi=300)
    plt.close()
    print(f"   -> Created {out_name}")
    return True

def run_automation():
    base_path = os.getcwd()
    print(f"--- Automated Raman Plotter (Style: Publication) ---")
    print(f"Scanning: {base_path}")
    
    fwhm, b_type = read_broadening_settings()
    print(f"Broadening: {'Lorentzian' if b_type == 'l' else 'Gaussian'}, FWHM = {fwhm:g} cm-1")

    count = 0
    broadened_count = 0
    avg_re = re.compile(r"^Raman_intensity_polarization_averaged_.+eV$")
    raw_complex_re = re.compile(r"^Raman_intensity_complex_(?!broadening_).+eV$")
    for root, dirs, files in os.walk(base_path):
        for fname in files:
            if raw_complex_re.match(fname):
                full_path = os.path.join(root, fname)
                try:
                    if write_broadened_spectrum(full_path, fwhm, b_type) is not None:
                        broadened_count += 1
                except Exception as error:
                    print(f"      Error broadening {fname}: {error}")
                print(f"Processing: {os.path.join(os.path.basename(root), fname)}")
                if process_and_plot(full_path, fwhm, b_type):
                    count += 1
            if avg_re.match(fname):
                full_path = os.path.join(root, fname)
                print(f"Processing: {os.path.join(os.path.basename(root), fname)}")
                if process_and_plot(full_path, fwhm, b_type):
                    count += 1

    print(f"\nSuccess! Generated {count} plots and {broadened_count} numerical broadened spectra.")

if __name__ == "__main__":
    run_automation()
