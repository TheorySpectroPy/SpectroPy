# SpectroPy

A Python toolkit for simulating Raman spectra from first-principles calculations. This package automates the setup of the finite displacement method for VASP, post-processes the results to calculate frequency-dependent Raman tensors, and generates plottable spectra for various experimental geometries.

## Key Features

- **Symmetry Reduction and Reconstruction**: Computes a symmetry-reduced set of atomic displacements and reconstructs the remaining dielectric derivatives at equivalent directions and atoms.
- **Symmetry Denoising**: Enforces atomic symmetry and reciprocity while retaining raw results.
- **Automated Setup**: Generates calculation directories and displaced `POSCAR` files.
- **Data Harvesting**: Automatically collects results from multiple `vasprun.xml` files.
- **Frequency-Dependent Tensors**: Calculates complex Raman tensors for user-specified laser energies.
- **Spectrum Simulation**: Produces plottable spectra for specific polarization geometries or summed polarization channels.
- **Mode Visualization**: Generates output files for visualizing phonon modes in **VMD** or **VESTA**.

## Requirements and Installation

Before you begin, ensure you have the following installed and available:

| Requirement | Purpose |
| --- | --- |
| **VASP** | Run displaced-structure optical calculations. A licensed installation is required; it is not included. |
| **Python 3.9+**, NumPy, PyYAML | Run SpectroPy; Python dependencies are installed automatically. |
| **Phonopy** | Supply Gamma-point phonons with eigenvectors, and symmetry information for reduced displacements. |
| **Matplotlib** | Plot broadened spectra. |
| **VMD / VESTA** (optional) | Visualize phonon modes. |

Install from a checkout:

```bash
python -m pip install .          # Core analysis
python -m pip install '.[all]'   # Also install plotting and Phonopy dependencies
spectropy --help
```

The dielectric reader accepts complete `vasprun.xml` files with either nine
components or the standard six-component dielectric format. All printed entries
are preserved independently. If transpose components are missing, SpectroPy
warns and fills them by reciprocity (for example, $\epsilon_{yx}=\epsilon_{xy}$,
without complex conjugation). These entries are inferred, not measured; this
cannot recover a nonreciprocal response. Malformed XML, missing diagonal entries,
or missing entire off-diagonal pairs still fail. There is no OUTCAR fallback.

## Workflow & Usage Guide

Follow these steps to calculate a Raman spectrum from a relaxed crystal structure.

```text
Relaxed structure + Gamma-point phonons
                 |
        spectropy displacements
                 |
       Run VASP in each directory
                 |
         spectropy derivatives
                 |
          spectropy spectrum
                 |
            spectropy plot
```

### Step 1: Prepare Your Directory

Create a new directory for your Raman calculation. Inside this directory, you will need:

| File | Contents |
| --- | --- |
| `CONTCAR` | Relaxed VASP structure. |
| `band.yaml` | Phonopy Gamma-point frequencies and mass-weighted eigenvectors. |
| `irreps.yaml` | Phonon symmetry labels for the plots. |
| `symmetry` | Phonopy YAML symmetry output. |
| `INCAR`, `KPOINTS`, `POTCAR` | VASP inputs to use for the finite displacement calculations. |
| `input.yaml` | Experimental geometry and SpectroPy settings, as below. |

Keep the calculation settings together in a file named `input.yaml`:

```yaml
incident_polarization: [1, 0, 0]
scattered_polarization: [0, 1, 0]
surface_normal: z
laser_energies: [0.00, 1.96] # eV: static response and one laser energy
symmetry: on # reduce atoms and directions using symmetry
displacement_basis: cartesian # basis for the candidate displacement vectors
displacement_amplitude: 0.01 # physical displacement length in Angstrom
difference_scheme: central # use positive/negative displacement pairs
site_symmetry_projection: true # denoise within the symmetry-allowed subspace
polarization: specific # use the polarization vectors above
broadening_hwhm: 1.0 # peak half-width in cm^-1
broadening_type: lorentzian # spectral line shape
```

Here, the geometry selects xy polarization with surface normal z. The energies
select static and 1.96 eV responses; `symmetry: on`, `cartesian`, and `0.01` select
symmetry-reduced directions with a 0.01 Å step. Central differences use ± pairs,
projection enables denoising, and Lorentzian peaks have a 1 cm⁻¹ half-width.
The step is illustrative, not convergence-tested. Each stage reads `input.yaml`.

### Step 2: Generate Displacements

Choose one `symmetry` entry in `input.yaml`. For this workflow, keep
`on` enabled; switch to `off` when debugging or comparing unreduced results.

```yaml
# symmetry: off    # Debugging: all atoms, three basis axes
symmetry: on  # Demonstration: symmetry-reduced atoms and directions
```

Then generate the structures:

```bash
spectropy displacements
```

**Input**: `CONTCAR`, `input.yaml`, and `symmetry` for the reduced example above.

The main generated files have this layout (`...` represents other displacements):

```text
calculation/
├── displacements.dat
├── ref_poscar.vasp
└── displacements/
    ├── ra_pos_atom1a/
    │   ├── POSCAR
    │   ├── INCAR   -> ../../INCAR
    │   ├── KPOINTS -> ../../KPOINTS
    │   └── POTCAR  -> ../../POTCAR
    └── ...
```

| File | Format and contents |
| --- | --- |
| `displacements.dat` | Plain text: selection metadata, atom indices, and fractional displacement vectors. |
| `ref_poscar.vasp`, `POSCAR` | VASP structure format: equilibrium reference or individual displaced structures. |
| `INCAR`, `KPOINTS`, `POTCAR` | Symlinks to your shared optical inputs, not copies. |

Displacement vectors in the text files are fractional even when candidate
directions are selected in the Cartesian basis; step lengths are in Å.

### Step 3: Run VASP

This is the main computational step. You must now run VASP in each of the
`displacements/ra_pos_*` subdirectories that were just created. This is your
responsibility and will depend on your specific computer or cluster environment.
Calculate the dielectric response without relaxing the displaced structures.

### Step 4: Calculate Dielectric Derivatives

Once all VASP jobs are complete:

```bash
spectropy derivatives
```

SpectroPy collects the displaced `vasprun.xml` files into `displacements/vasprun/`
and parses each once for all requested laser energies.

| Output (`<E>` is the laser energy in eV) | Format and contents |
| --- | --- |
| `epsilon_derivative_<E>` | Text: reconstructed, denoised dielectric derivatives for every atom, including symmetry-equivalent atoms and directions; three Cartesian 3×3 slices per atom, with real and imaginary parts. |
| `corrected_epsilon_derivative_<E>` | Text: explicitly named copy of the same denoised derivatives. |
| `raw_directional_derivative_<E>.json` | JSON: raw measured-direction derivatives, displacement vectors, step lengths, and sign provenance. |
| `dielectric_tensor_<E>` | Text: dielectric tensors sampled from the displacement calculations. |

For this example, the reconstructed files are `epsilon_derivative_0.00` and
`epsilon_derivative_1.96`, alongside their `corrected_` copies used by the
spectrum stage.

### Step 5: Calculate the Final Spectrum

```bash
spectropy spectrum
spectropy plot
```

SpectroPy combines the atomic derivatives with the phonon eigenvectors
(`band.yaml`) to calculate the mode Raman tensors and corresponding intensities.

| Output (`<E>` is the laser energy in eV) | Contents |
| --- | --- |
| `Raman_tensor_<E>eV` | Complex Raman tensor for each vibrational mode. |
| `Raman_intensity_complex_<E>eV` | Frequencies, intensities, and available symmetry labels for the chosen polarization. |
| `Raman_intensity_complex_broadening_<E>eV` | Unnormalized, two-column broadened spectrum. |

This example also saves `raw_` and `corrected_` tensor/intensity files.
The unprefixed files used by the plotter contain the denoised result.

![Example Raman Spectrum Plot](docs/images/Raman_plot_MoS2.png)

For the full settings reference, reconstruction equations, and other options,
see [Advanced Usage](docs.md).

## Additional Utilities

### Visualizing Phonon Modes

This script generates files to visualize the atomic motion of each phonon mode in either **VMD** or **VESTA**.

For VESTA, open `CONTCAR` and save a project file (e.g. `template.vesta`) in
your working directory. From that directory, run the installed helper:

```bash
python -m visualize_modes
```

| Prompt | Choice |
| --- | --- |
| Arrow scaling factor | For example, `0.5`; a negative value reverses the arrows. |
| Output format | `e` for VESTA, `v` for VMD, or `b` for both. |
| Template file | Your `.vesta` project, when using VESTA. |
| Supercell visualization | Choose whether to repeat the cell and, if so, enter dimensions such as `2 2 1`. |

The script will generate new folders (`VESTA_MODES` and/or `VMD_MODES`)
containing the visualization files for each mode. Open `.vesta` files in
VESTA; for VMD, go to `Extensions > Tk Console` and type:

```tcl
source VMD_MODES/mode_001.vmd
```

## Citation

If you use **SpectroPy** in your research, please cite our paper:

> **Raman Digital Twin of Monolayer Janus Transition Metal Dichalcogenides**
>
> *Johnathan Kowalski and Liangbo Liang*
>
> *ACS Applied Materials & Interfaces*, **2025**, *Article ASAP*
>
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
  url = {https://doi.org/10.1021/acsami.5c20316}
}
```

## License

SpectroPy is released under the [MIT License](LICENSE).
