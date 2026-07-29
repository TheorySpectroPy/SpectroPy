# SpectroPy

SpectroPy is a Python toolkit for predicting Raman spectra from first-principles
calculations with the Placzek approximation. It prepares finite-displacement
structures, uses Phonopy symmetry information to reduce and reconstruct the
required dielectric-tensor derivatives, and combines those derivatives with
phonon eigenvectors to produce Raman tensors, intensities, and broadened plots.


![Example Raman spectrum plot](docs/images/Raman_plot_MoS2.png)

## Features

- Generates finite-displacement `POSCAR` files from a relaxed `CONTCAR`.
- Optionally uses Phonopy to generate a symmetry-reduced displacement set.
- Prepares one calculation directory per displacement without prescribing a
  VASP execution method.
- Collects VASP `vasprun.xml` outputs and calculates complex,
  frequency-dependent dielectric-tensor derivatives.
- Reconstructs unmeasured directions and symmetry-equivalent atoms from
  phonon symmetry operations.
- Builds per-mode Raman tensors and intensities for specified and
  orientation-averaged experimental geometries.
- Produces Lorentzian- or Gaussian-broadened Raman plots.
- Generates VMD or VESTA phonon-mode visualizations.

## Installation

Install from a checkout:

```bash
python -m pip install .
```

Core dependencies are Python 3.9+, NumPy, and PyYAML. To use Phonopy-based
minimal displacements and plotting, install the optional dependencies:

```bash
python -m pip install '.[all]'
```

VASP is an external, proprietary prerequisite. You also need a Phonopy
calculation containing Gamma-point modes (`q-position: [0, 0, 0]`).

### `input` file

Each calculation directory may contain an `input` file. Its first three lines
retain the established MoS₂-style geometry format (incident polarization,
scattered polarization, then surface normal). Add `laser_energies:` to request
one or more dielectric-derivative energies. If the line is omitted, SpectroPy
uses `0.00 eV`, the static-limit calculation.

```text
1.0 0.0 0.0
1.0 0.0 0.0
z
laser_energies: 0.00 1.96 2.33
broadening_fwhm: 1.0
broadening_type: lorentzian
```

`broadening_fwhm` is in cm⁻¹ and defaults to `1.0`; `broadening_type`
defaults to `lorentzian` and may instead be `gaussian`.

## Command-line interface

The installed `spectropy` command runs existing workflow stages in the current
directory. Use `--workdir` (or `-C`) to run a stage elsewhere.

```bash
spectropy --help
spectropy displacements --help
spectropy displacements --mode minimal --workdir /path/to/calculation
```

| Command | Purpose |
| --- | --- |
| `spectropy displacements --mode full` | Create the full `+/-x`, `+/-y`, `+/-z` set for every atom and VASP-ready directories. |
| `spectropy displacements --mode atoms` | Create the full `+/-x`, `+/-y`, `+/-z` set for symmetry-inequivalent atoms only. |
| `spectropy displacements --mode minimal` | Create Phonopy's minimum symmetry-reduced displacement set. |
| `spectropy derivatives` | Process Phonopy symmetry, then evaluate dielectric derivatives at every laser energy in `input`. |
| `spectropy spectrum` | Combine `band.yaml` and derivatives into Raman tensors and intensities. |
| `spectropy plot` | Interactively broaden and plot generated Raman intensities. |

The CLI uses the same input/output conventions as the Python API. The optional
phonon-mode visualization helper remains available directly:

```bash
python src/visualize_modes.py
```

## Raman workflow

The following workflow starts from a relaxed crystal structure. The
[`examples/MoS2_Tutorial`](examples/MoS2_Tutorial/README.md) directory provides an example.

### 1. Prepare a calculation directory

Create a directory containing:

1. A relaxed VASP structure named `CONTCAR`.
2. Phonopy's Gamma-point `band.yaml` and, optionally, `irreps.yaml`.
3. Phonopy's YAML `symmetry` output when using symmetry-reduced displacements
   or derivative reconstruction.
4. The DFT inputs to reuse for every displacement: normally `INCAR`,
   `KPOINTS`, and `POTCAR` for VASP.

### 2. Generate displacement structures

```bash
# Full finite-difference set: every atom in +/- Cartesian directions.
spectropy displacements --mode full

# Full +/- Cartesian set, but only one atom from each symmetry-equivalent set.
# Requires Phonopy's `symmetry` file.
spectropy displacements --mode atoms

# Symmetry-reduced set (requires Phonopy).
spectropy displacements --mode minimal
```

Every mode writes `ref_poscar.vasp`, `displacements.dat`, and the
pipeline-compatible `atomic_displacement`, then creates `pos_atom*` and
`ra_pos_atom*/POSCAR` outputs. SpectroPy does not link or create DFT input
files; add the inputs required by your own DFT code to each directory. The
`atoms` mode requires Phonopy's YAML `symmetry` file and reconstructs the
remaining atoms during the derivative stage.

### 3. Run DFT calculations

Run your DFT code in every `ra_pos_atom*` directory. SpectroPy does not
choose your functional, k-point mesh, dielectric settings, executable, or job
submission strategy. For the VASP frequency-dependent workflow, each completed
directory must contain `vasprun.xml` with dielectric-function data at the
requested laser energy.

### 4. Process symmetry and dielectric derivatives

```bash
spectropy derivatives
```

This command processes Phonopy's `symmetry` file, then evaluates every
`laser_energies` value in `input`. It copies
`ra_pos_atom*/vasprun.xml` into `vasprun/`, then writes:

- `dielectric_tensor_<frequency>`: dielectric tensor for every displacement.
- `epsilon_derivative_<frequency>`: complex dielectric-tensor derivatives
  for all atoms, including symmetry-reconstructed data where applicable.

When `laser_energies` is absent, SpectroPy evaluates the default `0.00 eV`
case. This replaces the separate static-derivative workflow.

### 5. Calculate the Raman spectrum

```bash
spectropy spectrum
```

This command reads `band.yaml`, `irreps.yaml` when present, and the derivative
file for every `laser_energies` value in `input`. It prompts for the
experimental geometry unless `input` is provided, then writes:

- `Raman_tensor`: complex Raman tensor for each vibrational mode.
- `Raman_intensity_complex_<energy>eV`: frequency and intensity for the
  requested polarization geometry.
- `Raman_intensity_polarization_averaged_<energy>eV`: orientation-averaged
  (powder) intensities.

To predefine the geometry and laser energies, create `input` as follows:

```text
1.0 0.0 0.0   ! incident-light polarization
0.0 1.0 0.0   ! scattered-light polarization
z             ! surface-normal direction
laser_energies: 0.00 1.96 2.33 ! eV; 0.00 is the static-limit calculation
broadening_fwhm: 1.0 ! cm-1
broadening_type: lorentzian
```

### 6. Broaden and plot results

```bash
spectropy plot
```

The plotter reads `broadening_fwhm` (default `1.0 cm⁻¹`) and
`broadening_type` (`lorentzian` by default, or `gaussian`) from `input`. It
scans the current directory tree for Raman intensity files and writes
publication-style PNG plots alongside them. It also writes unnormalized,
two-column numerical spectra named
`Raman_intensity_complex_broadening_<energy>eV`, matching the Raman-workflow
file convention. If LaTeX is installed it is used for rendering; otherwise
Matplotlib's standard fonts are used.

## Phonon-mode visualization

`visualize_modes.py` generates VMD, VESTA, or both kinds of visualization
files from the phonon modes. For VESTA, first open `CONTCAR` in VESTA and save
a project file such as `template.vesta` in the calculation directory.

```bash
python visualize_modes.py
```

The script prompts for arrow scaling, output format, and—when needed—the VESTA
template. It creates `VESTA_MODES` and/or `VMD_MODES`. To load a VMD mode:

```tcl
source VMD_MODES/mode_001.vmd
```

## Citation

If you use SpectroPy in your research, please cite:

> **Raman Digital Twin of Monolayer Janus Transition Metal Dichalcogenides**
> Johnathan Kowalski and Liangbo Liang
> *ACS Applied Materials & Interfaces* (2025), Article ASAP.
> DOI: [10.1021/acsami.5c20316](https://doi.org/10.1021/acsami.5c20316)

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
