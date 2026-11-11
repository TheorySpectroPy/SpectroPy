# SpectroPy

A Python toolkit for simulating Raman spectra from first-principles calculations. This package automates the setup of the finite displacement method for VASP, post-processes the results to calculate frequency-dependent Raman tensors, and generates plottable spectra for various experimental geometries.

## Key Features

- **Automated Setup**: Generates calculation directories and displaced `POSCAR` files, with optional symmetry reduction.
- **Data Harvesting**: Automatically collects results from multiple `vasprun.xml` files.
- **Frequency-Dependent Tensors**: Calculates complex Raman tensors for user-specified laser energies.
- **Spectrum Simulation**: Produces plottable spectra for specific polarization geometries or summed polarization channels.
- **Symmetry Denoising**: Optionally enforces atomic symmetry and reciprocity while retaining raw results.
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

The dielectric reader requires a complete `vasprun.xml` containing all **nine**
named tensor components in its density-density dielectric block. Verify your
VASP build's output: version alone does not guarantee this format. Missing
transpose components are rejected, not inferred; there is no OUTCAR fallback.

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
| `irreps.yaml` (optional) | Phonon symmetry labels. |
| `symmetry` | Phonopy YAML symmetry output; required for reduced displacements or denoising. |
| `INCAR`, `KPOINTS`, `POTCAR` | VASP inputs to use for the finite displacement calculations. |
| `input` | Experimental geometry and SpectroPy settings, as below. |

To predefine the experimental geometry, create a file named `input`:

```text
1.0 0.0 0.0   ! polarization of the incident light
0.0 1.0 0.0   ! polarization of the scattered light
z             ! the surface normal (or out-of-plane) direction
laser_energies: 0.00 1.96
displacement_amplitude: 0.01
difference_scheme: central
site_symmetry_projection: false
broadening_hwhm: 1.0
```

### Step 2: Generate Displacements

```bash
spectropy displacements --mode full --basis cartesian
```

**Input**: `CONTCAR`, `input`

**Output**: `ref_poscar.vasp`, `displacements.dat`, `atomic_displacement`,
and `displacements/ra_pos_atom*/POSCAR`. An undisplaced copy is also written
to `displacements/ra_pos_equilibrium/POSCAR`.

Use `--mode minimal` for site-symmetry reduction. Details are below.

### Step 3: Prepare VASP Calculation Directories

Place your `INCAR`, `KPOINTS`, and `POTCAR` in every `displacements/ra_pos_*`
directory. SpectroPy prepares the structures, but does not populate these
electronic-structure inputs or select convergence settings.

### Step 4: Run VASP

This is the main computational step. You must now run VASP in each of the
`displacements/ra_pos_*` subdirectories that were just created. This is your
responsibility and will depend on your specific computer or cluster environment.
Calculate the dielectric response without relaxing the displaced structures.

### Step 5: Calculate Dielectric Derivatives

Once all VASP jobs are complete:

```bash
spectropy derivatives
```

SpectroPy collects the displaced `vasprun.xml` files into `displacements/vasprun/`
and parses each once for all requested laser energies.

### Step 6: Calculate the Final Spectrum

```bash
spectropy spectrum
spectropy plot
```

SpectroPy combines the atomic derivatives with the phonon eigenvectors
(`band.yaml`) to calculate the mode Raman tensors and corresponding intensities.

| Output (`<E>` is the laser energy in eV) | Contents |
| --- | --- |
| `dielectric_tensor_<E>` | Dielectric tensors for the displacements. |
| `raw_directional_derivative_<E>.json` | Raw derivatives in measured displacement directions, with direction vectors and sign provenance. |
| `epsilon_derivative_<E>` | Reconstructed Cartesian atomic derivatives; denoised when enabled. |
| `corrected_epsilon_derivative_<E>` | Additional denoised derivative file, when enabled. |
| `Raman_tensor_<E>eV` | Complex Raman tensor for each vibrational mode. |
| `Raman_intensity_complex_<E>eV` | Frequencies, intensities, and available symmetry labels for the chosen polarization. |
| `Raman_intensity_complex_broadening_<E>eV` | Unnormalized, two-column broadened spectrum. |

The atomic spectrum stage also writes `raw_` tensor/intensity files and, with
denoising enabled, `corrected_` files. Unprefixed files select the active result
used by the plotter. Existing results are not deleted when denoising is disabled.

![Example Raman Spectrum Plot](docs/images/Raman_plot_MoS2.png)

| Example | Starting point |
| --- | --- |
| [MoS2_Tutorial](examples/MoS2_Tutorial/README.md) | Bundled derivative data: calculate tensors and spectra. |
| [MoS2_From_Displacements](examples/MoS2_From_Displacements/README.md) | Detailed earlier-stage walkthrough; structure, phonon, and raw dielectric bundles are placeholders, not yet supplied. |
| [MoS2_Mode_Displacements](examples/MoS2_Mode_Displacements/README.md) | Displace collective Raman-active modes instead of individual atoms. |

## Advanced Usage

### Command-line reference

All subcommands accept `--workdir PATH` (or `-C PATH`, default `.`) and
`--help`. The top-level command also accepts `--version`.

| Command | Options and defaults |
| --- | --- |
| `spectropy displacements` | `--mode full` (`full`, `atoms`, `minimal`); `--basis cartesian` (`cartesian`, `fractional`). |
| `spectropy symmetry` | Write equivalent-atom mappings to `symmetry_operation_matrices`. |
| `spectropy derivatives` | `--difference central\|forward` overrides `input`; `--equilibrium-file PATH` overrides the default equilibrium XML; `--saved-basis` reuses validated reconstruction metadata. |
| `spectropy spectrum` | Read geometry and energies from `input`; prompt for geometry if no `input` file exists. |
| `spectropy plot` | Broaden and plot intensity files under the working directory. |
| `spectropy mode-displacements` | `--q 0.01`; `--modes N ...`; `--minimum-frequency 1.0`; `--output mode_displacements`. |
| `spectropy mode-tensors` | `--directory mode_displacements`; `--output mode_results`; `--source vasprun.xml`; `--energies E ...` (default: `input`). |

`--equilibrium-vasprun` is an alias for `--equilibrium-file`. Without either,
forward differences use `displacements/ra_pos_equilibrium/vasprun.xml`;
central differences do not use its response.

### Pre-defining Experimental Geometry and Settings

The first three lines of `input` contain the incident and scattered polarization
vectors and surface-normal axis (`x`, `y`, or `z`). Subsequent settings are:

| Setting | Default | Meaning |
| --- | --- | --- |
| `laser_energies` | `0.00` | One or more photon energies in eV; zero is the static limit. |
| `displacement_amplitude` | `0.03` | Physical displacement length in Å, in either basis. |
| `difference_scheme` | `central` | `central` or `forward`; set before generating displacements. |
| `dielectric_source` | `vasprun.xml` | Only supported dielectric source. |
| `site_symmetry_projection` | `false` | `true` enables atomic symmetry and reciprocity denoising. |
| `polarization` | `specific` | `specific` or `average`; the latter sums polarization channels, not a general powder orientational integral. |
| `cutoff` | Unset | When set, apply the 298 K Bose/phonon-frequency factor to modes at or above this threshold (cm⁻¹). Otherwise use reduced intensities. |
| `broadening_hwhm` | `1.0` | Half-width at half-maximum in cm⁻¹; FWHM is twice this value. |
| `broadening_type` | `lorentzian` | `lorentzian` or `gaussian` (also `l` or `g`). |

`full` and `atoms` accept three displacement lengths, one per selected basis
axis; `minimal` requires a common length. Mode and basis are CLI options, not
`input` settings. Accepted aliases are `laser_energy` and `frequency_cutoff`.
Legacy `broadening_fwhm` is converted to HWHM; do not specify both width keys.

### Automating Energy Scans

Request multiple energies in `input`, then run the analysis once:

```text
laser_energies: 0.00 0.10 0.20 0.30 1.50 1.96 2.33
```

```bash
spectropy derivatives
spectropy spectrum
spectropy plot
```

The nearest tabulated energy is selected without interpolation. Filenames use
two decimal places; choose energies that remain distinct at that precision.

### Displacement Selection and Derivatives

| Mode | Atoms | Directions |
| --- | --- | --- |
| `full` | Every atom | Three selected-basis axes. |
| `atoms` | One representative per atomic orbit | Three selected-basis axes. |
| `minimal` | One representative per atomic orbit | Site-symmetry-reduced axis or mixed directions. |

All three support Cartesian or fractional candidates. SpectroPy generates
displacements itself; Phonopy supplies symmetry information, not the
displacement list. Both bases use the same ordered candidate/orbit search.
Counts and finite-step errors can differ; the search does not guarantee
globally optimal conditioning. Minimal mode records its basis and transformations
in `fractional_basis.json` for either choice of candidate basis.

For a Cartesian unit direction $\mathbf u$, path $\delta\mathbf r(\mu)=s_\mu\mathbf u$,
and physical step $\Delta$, evaluated at equilibrium:

$$
\frac{\partial\epsilon}{\partial s_\mu}
=\sum_i u_i\frac{\partial\epsilon}{\partial r_i(\mu)},
\qquad
\frac{\partial\epsilon}{\partial s_\mu}\approx
\begin{cases}
\dfrac{\epsilon(\mu,+\Delta\mathbf u)-\epsilon(\mu,-\Delta\mathbf u)}{2\Delta}
& \text{central},\\[6pt]
\dfrac{\epsilon(\mu,+\Delta\mathbf u)-\epsilon(0)}{\Delta}
& \text{forward}.
\end{cases}
$$

All dielectric matrix components are Cartesian, even when the displacement
direction was selected in a fractional basis. Raw JSON records preserve these
directional responses, not a manually selected subset of Cartesian slices.
In minimal mode, symmetry may supply the negative sign; the JSON records this.

### Symmetry Denoising and Raman Tensors

For `atoms` and `minimal`, denoising fits the measured directional responses
within each representative's site-allowed, reciprocal tensor subspace before
transfer to equivalent atoms. For `full`, it projects the Cartesian field over
the complete space group, including equivalent-atom relationships. The latter
uses an equal-weight Cartesian norm, not an observation-space fit for
nonorthogonal fractional directions.

Reciprocity means $\partial\epsilon_{\alpha\beta}/\partial r_i
=\partial\epsilon_{\beta\alpha}/\partial r_i$, without complex conjugation;
do not impose it on general magneto-optical responses. No acoustic sum rule
is imposed, and symmetry-allowed numerical errors can remain.

Following the legacy derivative-file convention, stored atomic responses
include $V_{\rm cell}/(4\pi)$ and have units of Å². The mode contraction is

$$
R_{\alpha\beta}(k)=
\sum_{\mu,i}
\left[\frac{V_{\rm cell}}{4\pi}
\frac{\partial\epsilon_{\alpha\beta}}{\partial r_i(\mu)}\right]
\frac{e_i^k(\mu)}{\sqrt{m_\mu}},
\qquad
I_k^{\rm reduced}=
\left|\mathbf g_s^{\mathsf T}R(k)\mathbf g_i\right|^2.
$$

The volume factor is already stored and must not be applied twice. The
eigenvectors are mass-weighted and normalized. These are Placzek-response
derivatives, not a full resonant scattering calculation.

### Collective Mode Displacements

As an alternative to atomic displacements, use `CONTCAR`, Gamma-point
`band.yaml`, `irreps.yaml`, and `symmetry` to select Raman-active modes:

```bash
spectropy mode-displacements --q 0.01
# Run VASP externally in mode_displacements/ra_POSCAR_mode*/.
spectropy mode-tensors
```

Each mode, including each degenerate partner, receives its own explicit
positive and negative displacement. Automatic selection includes Raman-active
modes across polarizations with frequencies above `--minimum-frequency`
(default 1 cm⁻¹). `--modes 4 5` overrides activity and frequency filtering.

$$
\delta r_i(\mu)=\pm Q_k\frac{e_i^k(\mu)}{\sqrt{m_\mu}},
\qquad
R(k)\approx\frac{V_{\rm cell}}{4\pi}
\frac{\epsilon(+Q_k)-\epsilon(-Q_k)}{2Q_k}.
$$

`--q` is in Å$\sqrt{\mathrm{amu}}$, not a maximum atomic displacement.
The generated directory contains `mode_manifest.json`, `mode_index`,
`POSCAR_mode*`, and `ra_POSCAR_mode*/POSCAR`. Its `runRA` script stages available
VASP inputs; it does not launch calculations.

`mode-tensors` validates displaced geometries and requires both signs, central
differences, and `polarization: specific`. It writes
`raw_mode_tensors_<E>eV.npz`, `Raman_tensor_<E>eV`, and
`Raman_intensity_complex_<E>eV`. Each mode tensor is extracted directly from
its own pair; atomic denoising is a separate workflow.

### Publication-Quality Plotting

`spectropy plot` generates Raman spectra plots with Lorentzian or Gaussian
broadening and labels prominent peaks with their available mode symmetries.
It scans the working directory tree, using the selected unprefixed intensities.
PNG plots are normalized for display; numerical broadened spectra are not.
LaTeX is used when installed, with a Matplotlib-font fallback.

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
