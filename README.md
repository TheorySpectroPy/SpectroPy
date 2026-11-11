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
| `irreps.yaml` (optional) | Phonon symmetry labels. |
| `symmetry` | Phonopy YAML symmetry output; required for reduced displacements or denoising. |
| `INCAR`, `KPOINTS`, `POTCAR` | VASP inputs to use for the finite displacement calculations. |
| `input` | Experimental geometry and SpectroPy settings, as below. |

Keep the calculation settings together in a file named `input`:

```text
incident_polarization: 1 0 0
scattered_polarization: 0 1 0
surface_normal: z
laser_energies: 0.00 1.96        ! eV: static response and one laser energy
displacement_mode: minimal     ! reduce atoms and directions using symmetry
displacement_basis: cartesian  ! basis for the candidate displacement vectors
displacement_amplitude: 0.01   ! physical displacement length in Angstrom
difference_scheme: central     ! use positive/negative displacement pairs
site_symmetry_projection: true ! denoise within the symmetry-allowed subspace
polarization: specific         ! use the polarization vectors above
broadening_hwhm: 1.0           ! peak half-width in cm^-1
broadening_type: lorentzian    ! spectral line shape
```

Here, the geometry selects xy polarization with surface normal z. The energies
select static and 1.96 eV responses; `minimal`, `cartesian`, and `0.01` select
symmetry-reduced directions with a 0.01 Å step. Central differences use ± pairs,
projection enables denoising, and Lorentzian peaks have a 1 cm⁻¹ half-width.
The step is illustrative, not convergence-tested. Each stage reads `input`;
omitted settings use the defaults below.

### Step 2: Generate Displacements

Choose one active `displacement_mode` entry in `input`. For this workflow, keep
`minimal` enabled; switch to `full` when debugging or comparing unreduced results.

```text
# displacement_mode: full    # Debugging: all atoms, three basis axes
displacement_mode: minimal  # Demonstration: symmetry-reduced atoms and directions
```

Then generate the structures:

```bash
spectropy displacements
```

**Input**: `CONTCAR`, `input`, and `symmetry` for the reduced example above.

| Output file | Format and contents |
| --- | --- |
| `ref_poscar.vasp` | VASP POSCAR format: equilibrium lattice, species, and fractional atomic positions. |
| `displacements.dat` | Plain text: mode/basis header, displacement count, then atom labels, 1-based indices, and fractional displacement vectors; ends with the displaced-atom list. |
| `atomic_displacement` | Fortran-compatible text layout: physical step lengths, displaced-atom list, and indexed fractional displacement vectors. |
| `fractional_basis.json` (minimal only) | JSON: reconstruction bases, inverse matrices, and symmetry mappings for each representative atom. |
| `displacements/pos_atom*` | VASP POSCAR format: one structure per displacement. |
| `displacements/ra_pos_atom*/POSCAR` | The same displaced structures, placed in individual VASP calculation directories. |
| `displacements/ra_pos_equilibrium/POSCAR` | Unchanged copy of the original `CONTCAR`, for the undisplaced optical calculation. |

Displacement vectors in the text files are fractional even when candidate
directions are selected in the Cartesian basis; step lengths are in Å.

Use `displacement_mode: full` to calculate every atom along three basis axes;
`atoms` reduces only the atoms, and `minimal` also reduces directions.

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

### Input Settings

Settings are named entries in `input`; their order does not matter:

| Setting | Default | Meaning |
| --- | --- | --- |
| `incident_polarization` | Required for spectra | Three Cartesian components of the incoming polarization vector. |
| `scattered_polarization` | Required for spectra | Three Cartesian components of the detected polarization vector. |
| `surface_normal` | Required for spectra | Surface-normal axis: `x`, `y`, or `z`. |
| `laser_energies` | `0.00` | One or more photon energies in eV; zero is the static limit. |
| `displacement_mode` | `full` | `full`, `atoms`, or `minimal`; see the selection table below. |
| `displacement_basis` | `cartesian` | `cartesian` or `fractional` candidate directions. |
| `displacement_amplitude` | `0.03` | Physical displacement length in Å, in either basis. |
| `difference_scheme` | `central` | `central` or `forward` derivative analysis. Generation retains its signed displacement set for either scheme. |
| `dielectric_source` | `vasprun.xml` | Only supported dielectric source. |
| `site_symmetry_projection` | `false` | `true` enables atomic symmetry and reciprocity denoising. |
| `polarization` | `specific` | `specific` or `average`; the latter sums polarization channels, not a general powder orientational integral. |
| `cutoff` | Unset | When set, apply the 298 K Bose/phonon-frequency factor to modes at or above this threshold (cm⁻¹). Otherwise use reduced intensities. |
| `broadening_hwhm` | `1.0` | Half-width at half-maximum in cm⁻¹; FWHM is twice this value. |
| `broadening_type` | `lorentzian` | `lorentzian` or `gaussian` (also `l` or `g`). |
| `equilibrium_file` | Automatic | Optional equilibrium XML path for forward differences; otherwise `displacements/ra_pos_equilibrium/vasprun.xml`. |
| `saved_basis` | `false` | Reuse validated `fractional_basis.json` when present instead of rebuilding the basis. |

`full` and `atoms` accept three displacement lengths, one per selected basis
axis; `minimal` requires a common length. Accepted aliases are `laser_energy` and `frequency_cutoff`.
Legacy `broadening_fwhm` is converted to HWHM; do not specify both width keys.

Paths are relative to the working directory. Changing `input` does not modify
existing structures or results: regenerate displacements in a fresh directory
to change their length, selection, or basis. Derivative analysis uses the
recorded displacement metadata, not a newly edited selection setting.

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

Add these settings to the same `input` file:

```text
mode_q: 0.01                         ! normal-coordinate amplitude, Angstrom sqrt(amu)
mode_indices: auto                   ! select Raman-active modes automatically
mode_minimum_frequency: 1.0          ! exclude modes at or below this cm^-1
mode_displacements_directory: mode_displacements
mode_results_directory: mode_results
```

These are the defaults. Set `mode_indices: 6 7`, for example, to select specific
1-based mode indices instead of applying the activity and frequency filters.
The directory settings control both generation and subsequent tensor extraction.

```bash
spectropy mode-displacements
# Run VASP externally in mode_displacements/ra_POSCAR_mode*/.
spectropy mode-tensors
```

Each mode, including each degenerate partner, receives its own explicit
positive and negative displacement. Automatic selection includes Raman-active
modes across polarizations with frequencies above `mode_minimum_frequency`
(default 1 cm⁻¹).

$$
\delta r_i(\mu)=\pm Q_k\frac{e_i^k(\mu)}{\sqrt{m_\mu}},
\qquad
R(k)\approx\frac{V_{\rm cell}}{4\pi}
\frac{\epsilon(+Q_k)-\epsilon(-Q_k)}{2Q_k}.
$$

`mode_q` is in Å$\sqrt{\mathrm{amu}}$, not a maximum atomic displacement;
it is independent of the atomic `displacement_amplitude` setting.
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

### Optional Command-line Overrides

The plain commands above are sufficient when settings are in `input`.
For one-off changes, an explicit CLI option overrides the corresponding setting
for that invocation only; it does not rewrite `input`.

**Precedence: explicit CLI option → `input` setting → built-in default.**

| Command | Optional overrides → `input` setting |
| --- | --- |
| `spectropy displacements` | `--mode` → `displacement_mode`; `--basis` → `displacement_basis`. |
| `spectropy derivatives` | `--difference` → `difference_scheme`; `--equilibrium-file` → `equilibrium_file`; `--saved-basis` / `--no-saved-basis` → `saved_basis`. |
| `spectropy mode-displacements` | `--q` → `mode_q`; `--modes` → `mode_indices`; `--minimum-frequency` → `mode_minimum_frequency`; `--output` → `mode_displacements_directory`. |
| `spectropy mode-tensors` | `--directory` → `mode_displacements_directory`; `--output` → `mode_results_directory`; `--source` → `dielectric_source`; `--energies` → `laser_energies`. |

All subcommands accept `--workdir PATH` (or `-C PATH`, default `.`) to choose
the directory containing `input`, and `--help` to list options. The top-level
command accepts `--version`. `--equilibrium-vasprun` is an alias for
`--equilibrium-file`. `spectropy symmetry` writes equivalent-atom mappings to
`symmetry_operation_matrices`; `spectrum` and `plot` use their `input` settings
without additional scientific CLI options.

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
