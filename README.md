# SpectroPy

A Python toolkit for simulating Raman spectra from first-principles calculations. This package automates the setup of the finite displacement method for VASP, post-processes the results to calculate frequency-dependent Raman tensors, and generates plottable spectra for various experimental geometries.

## Key Features

- **Automated Setup**: Generates all necessary directories and displaced `POSCAR` files for finite displacement calculations.
- **Data Harvesting**: Automatically collects results from multiple `vasprun.xml` files.
- **Frequency-Dependent Tensors**: Calculates complex Raman tensors for a user-specified laser frequency.
- **Spectrum Simulation**: Produces final, plottable spectra for both specific polarization geometries and orientation-averaged cases.

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
- Automatically collects all `vasprun.xml` files from the `raman_poscar_*` directories and copies them to a new `AXML` folder.  
- Calculates the derivative of the dielectric tensor with respect to each atomic displacement.  

**Outputs**:  
- `EPSILON_<freq>.dat`: A log file containing the dielectric tensor for each displacement.  
- `dielectric_derivatives_<freq>.dat`: The final calculated atomic Raman tensors for the specified frequency.

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
- `Raman_intensity_specific.dat`: A plottable (Frequency vs. Intensity) spectrum for your chosen polarization geometry.  
- `Raman_intensity_averaged.dat`: A plottable spectrum simulating a randomly oriented or powder sample.

---

## Advanced Usage

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

To avoid the interactive prompts in `calculate_spectrum.py`, you can create a file named `input` with the following format:

```
1.0 0.0 0.0   ! polarization of the incident light
0.0 1.0 0.0   ! polarization of the scattered light
z             ! the surface normal (or out-of-plane) direction
```


