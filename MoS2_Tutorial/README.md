# Tutorial: Monolayer MoS₂ Raman Spectrum

This folder contains a complete set of input files and pre-calculated intermediate data for simulating the Raman spectrum of **Monolayer MoS₂**.

You can use this tutorial to verify that your installation is working correctly and to familiarize yourself with the workflow.

## Included Files

We have provided the essential files, so you can skip the VASP steps:

### Inputs

- `CONTCAR`, `INCAR`, `KPOINTS`: Standard VASP inputs for the simulation.
- `band.yaml`: Phonon modes calculated by Phonopy.

### Pre-Calculated Data

- `dielectric_derivatives_1.96.dat`: The Raman tensors calculated from VASP (so you don't have to run the supercomputer jobs).

## Running the Tutorial

### Step 1: Generate Displacements

First, test the setup script to see how it handles the structure.

Run from this directory:

```bash
python ../../create_displacements.py
```

**Expected Output:**  
A new file `displacements.dat` will appear, listing the required atomic movements.

---

### Step 2: Calculate the Spectrum

Since we already provided the `dielectric_derivatives_1.96.dat` file, you can skip the VASP calculations and go straight to the physics analysis.

```bash
python ../../calculate_spectrum.py
```

**Prompts to answer:**

- **Laser Energy:** Enter `1.96` (to match the provided data file).
- **Geometry:** You can choose specific geometries (like z direction), or just press Enter for defaults.

**Expected Output:**  
The script will combine the provided derivatives with the `band.yaml` info to create:

- `Raman_intensity_specific.dat`
- `Raman_intensity_averaged.dat`

---

### Step 3: Visualize the Result

Now, generate the publication-quality plot.

```bash
python ../../generate_raman_plots.py
```

**Prompts to answer:**

- **FWHM:** Press Enter (default `5.0`).
- **Broadening:** Press Enter (default `Lorentzian`).

**Final Result:**  
A file named `Raman_plot_styled.png` will be generated. It should look identical to the reference image.

![Example Raman Spectrum Plot](../docs/images/Raman_plot_MoS2.png)
