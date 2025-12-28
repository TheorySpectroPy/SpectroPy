# SpectroPy

A Python toolkit for simulating Raman spectra from first-principles calculations. This package automates the setup of the finite displacement method for VASP, post-processes the results to calculate frequency-dependent Raman tensors, and generates plottable spectra for various experimental geometries.

## Key Features

- **Automated Setup**: Generates all necessary directories and displaced `POSCAR` files for finite displacement calculations.
- **Data Harvesting**: Automatically collects results from multiple `vasprun.xml` files.
- **Frequency-Dependent Tensors**: Calculates complex Raman tensors for a user-specified laser frequency.
- **Spectrum Simulation**: Produces final, plottable spectra for both specific polarization geometries and orientation-averaged cases.
- **Mode Visualization**: Generates output files for visualizing phonon modes in **VMD** or **VESTA**.

---

## Requirements

Before you begin, ensure you have the following installed and available:

- A working installation of **VASP**.
- **Python 3.x** with the following packages:
    ```bash
    pip install numpy pyyaml
    ```
- The results of a **Phonopy** calculation for your material at the Gamma point (`q-position: [0, 0, 0]`).

---

## Workflow & Usage Guide

Follow these steps to calculate a Raman spectrum from a relaxed crystal structure.

### Step 1: Prepare Your Directory

Create a new directory for your Raman calculation. Inside this directory, you will need:

1. A relaxed VASP structure file named `CONTCAR`.
2. The results of your Phonopy calculation: `band.yaml` and (optionally) `irreps.yaml`.
3. The three standard VASP input files you intend to use for the finite displacement calculations: `INCAR`, `KPOINTS`, and `POTCAR`.

### Step 2: Generate Displacements

Run the `create_displacements.py` script. This will read your `CONTCAR` and create a "recipe" of all required atomic movements.

```bash
python create_displacements.py
```

**Input**: `CONTCAR`  
**Output**: `displacements.dat` (A list of every atom to be displaced and the displacement vector).

### Step 3: Prepare VASP Calculation Directories

Run the `prepare_vasp_inputs.py` script. This script reads `displacements.dat` and automatically creates a dedicated subdirectory for each VASP calculation.

```bash
python prepare_vasp_inputs.py
```

**Input**: `CONTCAR`, `displacements.dat`  
**Action**:  
- Generates `poscar_<atom_id>` files for each displacement.  
- Creates a shell script `setup_vasp_calcs.sh`.  
- Automatically executes `setup_vasp_calcs.sh`, which creates all `raman_poscar_*` subdirectories and populates them with the correct `POSCAR` and symbolic links to your `INCAR`, `KPOINTS`, and `POTCAR`.

### Step 4: Run VASP

This is the main computational step. You must now run VASP in each of the `raman_poscar_*` subdirectories that were just created. This is your responsibility and will depend on your specific computer or cluster environment.

### Step 5: Calculate Dielectric Derivatives

Once all VASP jobs are complete, run the `calculate_dielectric_derivatives.py` script. It will automatically find and process all the `vasprun.xml` files.

```bash
python calculate_dielectric_derivatives.py
```

**Action**:  
- Prompts you to enter the laser frequency in eV.  
- Automatically collects all `vasprun.xml` files from the `raman_poscar_*` directories and copies them to a new `vasprun` folder.  
- Calculates the derivative of the dielectric tensor with respect to each atomic displacement.  

**Outputs**:  
- `EPSILON_<freq>.dat`: A log file containing the dielectric tensor for each displacement.  
- `dielectric_derivatives_<freq>.dat`: The final calculated atomic Raman tensors (dielectric tensor derivatives with respect to the atomic displacements) for the specified frequency.

### Step 6: Calculate the Final Spectrum

Finally, run the `calculate_spectrum.py` script to generate the plottable results.

```bash
python calculate_spectrum.py
```

**Action**:  
- Prompts you for the experimental geometry (polarization of incident/scattered light).  
- Combines the atomic Raman tensors (`dielectric_derivatives_<freq>.dat`) with the phonon eigenvectors (`band.yaml`).  
- Calculates the final Mode Raman Tensors and corresponding intensities.  

**Outputs**:  
- `Raman_tensor.dat`: The full complex Raman tensor for each vibrational mode.  
- `Raman_intensity_specific.dat`: Frequency, intensity, and irreducible representations for your chosen polarization geometry.  
- `Raman_intensity_averaged.dat`: Frequency, intensity, and irreducible representations for a randomly oriented or powder sample.

---

## Advanced Usage

### Publication-Quality Plotting (`generate_raman_plots.py`)

This script automates the generation of publication-quality Raman spectra plots (with Lorentzian/Gaussian broadening) for every calculation in your project directory.

![Example Raman Spectrum Plot](docs/images/Raman_plot_MoS2.png)

**Features:**
- **Batch Processing:** Automatically scans all subdirectories for `Raman_intensity_specific.dat`.
- **Publication Style:** Uses LaTeX formatting (so you will need it installed), inward ticks, and professional styling.
- **Smart Labeling:** Automatically labels prominent peaks with their mode symmetry (e.g., $E_{2g}$), using arrows to point to the curve.

**How to use:**
1. Place `generate_raman_plots.py` in your top-level project folder (e.g., `TMDs`).
2. Run the script:
   ```bash
   python generate_raman_plots.py

### Automating Energy Scans

To calculate the Raman response over a range of laser frequencies, you can use the provided `run_energy_scan.sh` script.

1. Make the script executable:
     ```bash
     chmod +x run_energy_scan.sh
     ```
2. Run the script. It will prompt you for a start energy, end energy, and increment.
     ```bash
     ./run_energy_scan.sh
     ```

This will automatically call `calculate_dielectric_derivatives.py` for each energy step, generating separate `dielectric_derivatives_<freq>.dat` and `EPSILON_<freq>.dat` files for each one. You can then run `calculate_spectrum.py` on each of these files.

### Pre-defining Experimental Geometry

To avoid some of the interactive prompts in `calculate_spectrum.py`, you can create a file named `input` with the following format:

```
1.0 0.0 0.0   ! polarization of the incident light
0.0 1.0 0.0   ! polarization of the scattered light
z             ! the surface normal (or out-of-plane) direction
```

# Additional Utilities

## Visualizing Phonon Modes (`visualize_modes.py`)

This script generates files to visualize the atomic motion of each phonon mode in either **VMD** or **VESTA**.

---

### How to use:

#### For VESTA (Recommended):

1. Open your `CONTCAR` file in **VESTA**.  
2. Save the file as a VESTA project file (e.g., `template.vesta`) in your working directory.

---

### Run the script:
```bash
python visualize_modes.py
```

### Answer the prompts:

- **Arrow scaling factor:** A number to control the arrow length (e.g., `0.5`). A negative value can be used to reverse the arrow direction.
- **Select output format:** Type `e` for **VESTA**, `v` for **VMD**, or `b` for **Both**.
- **Enter template file:** (If you chose VESTA) Type the name of the file you saved in step 1 (e.g., `template.vesta`).

---

### Output

The script will generate new folders (`VESTA_MODES` and/or `VMD_MODES`) containing the visualization files for each mode.

- **`.vesta` files:** Open these in **VESTA** to see the structure and arrows automatically.  
- **`.vmd` files:** Open **VMD**, go to `Extensions > Tk Console`, and type:
  ```tcl
  source VMD_MODES/mode_001.vmd
  ```
  to load the arrows for mode 1.


## Citation

If you use **SpectroPy** in your research, please cite our paper:

> **Raman Digital Twin of Monolayer Janus Transition Metal Dichalcogenides**
> *Johnathan Kowalski and Liangbo Liang*
> *ACS Applied Materials & Interfaces*, **2025**, *Article ASAP*
> **DOI:** [10.1021/acsami.5c20316](https://doi.org/10.1021/acsami.5c20316)

### BibTeX

```bibtex
@article{SpectroPy_Janus2025,
  author = {Kowalski, Johnathan and Liang, Liangbo},
  title = {Raman Digital Twin of Monolayer Janus Transition Metal Dichalcogenides},
  journal = {ACS Applied Materials \& Interfaces},
  year = {2025},
  note = {Article ASAP},
  doi = {10.1021/acsami.5c20316},
  url = {[https://doi.org/10.1021/acsami.5c20316](https://doi.org/10.1021/acsami.5c20316)}
}

