import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import sys
import re

# --- 1. Publication-Style Configuration ---
try:
    # Uses Helvetica/Sans-Serif
    rc('text', usetex=True)
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica', 'Arial']})
except:
    print("Notice: LaTeX not found. Using standard Matplotlib fonts.")
    rc('font', family='sans-serif')

# --- 2. Math & Helper Functions ---

def gaussian(x, center, amplitude, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return amplitude * np.exp(-((x - center)**2) / (2 * sigma**2))

def lorentzian(x, center, amplitude, fwhm):
    hwhm = fwhm / 2.0
    return amplitude * (hwhm**2) / ((x - center)**2 + hwhm**2)

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

# --- 3. The Plotting Core ---

def process_and_plot(input_file, fwhm=5.0, b_type='l'):
    try:
        # Load Data: Skip header, read columns 0, 1, 2
        raw_data = np.loadtxt(input_file, skiprows=1, dtype=str)
        if raw_data.size == 0: return
        if raw_data.ndim == 1: raw_data = raw_data.reshape(1, -1)
        
        freqs = raw_data[:, 0].astype(float)
        intensities = raw_data[:, 1].astype(float)
        modes = raw_data[:, 2]
    except Exception as e:
        print(f"      Error reading {os.path.basename(input_file)}: {e}")
        return

    # Generate Broadened Curve
    # Create x-axis with buffer
    x_dense = np.linspace(max(0, min(freqs)-50), max(freqs)+50, 2000)
    y_dense = np.zeros_like(x_dense)
    
    for f, i in zip(freqs, intensities):
        if b_type == 'g':
            y_dense += gaussian(x_dense, f, i, fwhm)
        else:
            y_dense += lorentzian(x_dense, f, i, fwhm)

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
    
    # Save Output
    out_name = os.path.join(os.path.dirname(input_file), f"Raman_plot_styled.png")
    plt.savefig(out_name, dpi=300)
    plt.close()
    print(f"   -> Created {out_name}")

# --- 4. Automation Logic ---

def run_automation():
    base_path = os.getcwd()
    print(f"--- Automated Raman Plotter (Style: Publication) ---")
    print(f"Scanning: {base_path}")
    
    # Inputs
    val_fwhm = input("Enter FWHM (cm-1) [default 5.0]: ")
    fwhm = float(val_fwhm) if val_fwhm else 5.0
    
    val_type = input("Broadening [L]orentzian or [G]aussian [default L]: ").lower()
    b_type = val_type if val_type in ['l', 'g'] else 'l'

    count = 0
    # Walk Directory
    for root, dirs, files in os.walk(base_path):
        target_file = "Raman_intensity_specific.dat"
        
        if target_file in files:
            full_path = os.path.join(root, target_file)
            print(f"Processing: {os.path.basename(root)}")
            process_and_plot(full_path, fwhm, b_type)
            count += 1

    print(f"\nSuccess! Generated {count} plots.")

if __name__ == "__main__":
    run_automation()