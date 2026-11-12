# SpectroPy Advanced Usage

See the [basic workflow](README.md#workflow--usage-guide) for a short introduction.

## Input Settings

Settings are grouped YAML mappings in `input.yaml`; their order does not matter.
The table uses dotted paths (for example, `symmetry.denoise` means `denoise`
inside `symmetry`), not literal dotted YAML keys. See the basic example.
Use lists for vectors and multiple values, `true`/`false` for booleans, and
`#` for comments. Duplicate and unknown keys are rejected.
`symmetry.reduction: on` and `symmetry.reduction: off` are string choices, not booleans; quotes are
optional in SpectroPy. The loader uses only `true` and `false` as boolean words.

| Setting | Default | Meaning |
| --- | --- | --- |
| `polarization.incident` | Required for spectra | Three Cartesian components of the incoming polarization vector. |
| `polarization.scattered` | Required for spectra | Three Cartesian components of the detected polarization vector. |
| `polarization.surface_normal` | Required for spectra | Surface-normal axis: `x`, `y`, or `z`. |
| `spectrum.energies` | `0.00` | A scalar, list, or inclusive `{start, stop, step}` energy range in eV; zero is the static limit. |
| `symmetry.reduction` | `off` | `on` (minimal), `off` (full), or `atoms` (atoms only); see the selection table below. |
| `displacements.basis` | `cartesian` | `cartesian` or `fractional` candidate directions. |
| `displacements.amplitude` | `0.03` | Physical displacement length in Å, in either basis. |
| `displacements.difference` | `central` | `central` or `forward` derivative analysis. Only `forward` generates an equilibrium calculation directory; the signed displacement set is unchanged. |
| `symmetry.denoise` | `false` | `true` enables atomic symmetry and reciprocity denoising. |
| `polarization.mode` | `specific` | `specific` or `average`; the latter sums polarization channels, not a general powder orientational integral. |
| `spectrum.cutoff` | Unset | Apply the inverse phonon-frequency factor only to positive modes at or above this threshold (cm⁻¹). Unset disables the denominator; modes below the cutoff retain a denominator factor of 1, not zero intensity. |
| `spectrum.temperature` | `298` | Temperature in kelvin for the independent Stokes Bose factor, applied to all positive-frequency modes regardless of cutoff. At 0 K this factor is 1. Zero/negative-frequency modes receive neither correction. |
| `broadening.hwhm` | `1.0` | Half-width at half-maximum in cm⁻¹; FWHM is twice this value. |
| `broadening.type` | `lorentzian` | `lorentzian` or `gaussian` (also `l` or `g`). |
| `displacements.equilibrium_file` | Automatic | Optional equilibrium XML path for forward differences; otherwise `displacements/ra_pos_equilibrium/vasprun.xml`. |
| `symmetry.saved_basis` | `false` | For `on`, generate `fractional_basis.json` and reuse it after validation; otherwise reconstruct without a cache. |

`off` and `atoms` accept a list of three displacement lengths, one per selected
basis axis, e.g. `displacements.amplitude: [0.01, 0.02, 0.03]`; `on` requires
a common length.

The former positional `input` file and flat `input.yaml` settings are no longer
accepted. Rewrite them using the grouped schema; renaming alone is insufficient.
No scientific inputs or completed calculations are migrated automatically.
The dielectric source is always `vasprun.xml`; there is no source-selection knob.

Paths are relative to the working directory. Changing `input.yaml` does not modify
existing structures or results: regenerate displacements in a fresh directory
to change their length, selection, or basis. Derivative analysis uses the
recorded displacement metadata, not a newly edited selection setting.

## Auxiliary Files and Optional Behavior

| File | Purpose |
| --- | --- |
| `fractional_basis.json` | Optional JSON reconstruction cache, generated and used only when both `symmetry.saved_basis: true` and `symmetry.reduction: on` are set. |
| `atomic_displacement` | Fortran-format export of the displacement list; generated, but not read by SpectroPy. |
| `displacements/pos_atom*` | Standalone POSCAR exports duplicating the calculation-directory structures; not read by SpectroPy. |

Only `displacements.dat` and the per-calculation `POSCAR` files are essential
generated inputs. Keep the original `CONTCAR` and symmetry/phonon inputs too.
The tensor stage registers phonon atoms against the reconstruction's `CONTCAR`.
`irreps.yaml` supplies labels and can be omitted for atomic-tensor calculations.

For forward differences, set `displacements.difference: forward` before generating
displacements. This also creates `displacements/ra_pos_equilibrium/POSCAR`,
copied from the reference structure, with the same input symlinks. Run VASP
there with matching optical settings; analysis automatically reads its
`vasprun.xml`, or an existing response specified by `displacements.equilibrium_file`.
Central differences do not create or update this directory. Existing
equilibrium calculations are left untouched, not deleted.

Input symlinks are created for supplied files. Missing inputs produce a warning;
matching links are retained, and conflicting files or links are never overwritten.
Supply any missing inputs before running VASP.

## Optional Per-atom and Debug Outputs

By default, `spectropy tensors` writes only `Raman_tensor_<E>eV` analysis results.
Atomic reconstruction, optional denoising, and phonon contraction happen in memory.
The existing XML collection cache in `displacements/vasprun/` is still used.

To inspect the final per-atom derivatives, enable:

```yaml
output:
  per_atom:
    enabled: true
    raw_directional_derivatives: false
```

| Setting | Default | Additional output |
| --- | --- | --- |
| `output.per_atom.enabled` | `false` | `epsilon_derivative_<E>`: final real/imaginary Cartesian slices for every atom; `dielectric_tensor_<E>`: the sampled dielectric responses. |
| `output.per_atom.raw_directional_derivatives` | `false` | `raw_directional_derivative_<E>`: labeled text containing measured directional responses for debugging. Requires `enabled: true`. |

The final per-atom file is denoised when `symmetry.denoise: true`; otherwise it
contains the unprojected reconstructed field. There is no duplicate
`corrected_epsilon_derivative` file. The raw text export preserves measured
directions, steps, sign provenance, and separate real/imaginary 3×3 matrices; it does not
pretend a mixed displacement is a directly measured Cartesian slice.
Like the other tensor text outputs, it uses full-precision numerical entries.
It is an inspection output, not an input to the spectrum stage. Older JSON
exports are not read or deleted automatically.

`dielectric_tensor_<E>` contains dimensionless Cartesian real/imaginary 3×3
matrices for each directly calculated signed displacement, with its atom and
signed Cartesian displacement vector in Å. Forward differences additionally
include the equilibrium response. These are the nearest-grid XML samples,
before denoising and without the derivative files' volume factor. All nine
reader components are written independently at full precision; when a source
XML prints only six, transpose completion still triggers the usual warning.
No additional XML reads are needed for this export.

`spectropy spectrum` reads only the final mode tensor files, so per-atom/debug
exports are never required. Changing denoising or energy settings requires rerunning
`spectropy tensors`; changing polarization requires only `spectrum` and `plot`.
Old raw/corrected exports are left untouched and are not selected automatically.
Tensor files retain full numerical precision between the two stages; no duplicate
unqualified `Raman_tensor` alias or raw/corrected spectrum copies are produced.

## Automating Energy Scans

For a regular scan, use an inclusive range in `input.yaml`:

```yaml
spectrum:
  energies: {start: 0.0, stop: 4.0, step: 0.1}
```

This selects 41 energies, including 0.0 and 4.0 eV. Both endpoints must lie
on the step grid. The step must be positive, and a range is limited to 100000
values. Expansion avoids cumulative floating-point step drift. For selected
energies instead, use a list:

```yaml
spectrum:
  energies: [0.00, 0.10, 0.20, 0.30, 1.50, 1.96, 2.33]
```

```bash
spectropy tensors
spectropy spectrum
spectropy plot
```

The nearest tabulated energy is selected without interpolation; all requested
energies must lie within the XML grid. Filenames use two decimal places, and
settings with distinct energies that would overwrite the same filename are
rejected. The settings writer stores expanded lists, which preserve the same
energy selection as a range.

## Displacement Selection and Derivatives

| `symmetry.reduction` | Atoms | Directions |
| --- | --- | --- |
| `off` | Every atom | Three selected-basis axes. |
| `atoms` | One representative per atomic orbit | Three selected-basis axes. |
| `on` | One representative per atomic orbit | Site-symmetry-reduced axis or mixed directions. |

`symmetry.reduction` controls displacement selection, independently of denoising:
`symmetry.reduction: off` still allows `symmetry.denoise: true` when comparing
the full set against a symmetry-denoised result. The existing algorithm names
in `displacements.dat` metadata are retained so completed datasets remain readable.

All three support Cartesian or fractional candidates. SpectroPy generates
displacements itself; Phonopy supplies symmetry information, not the
displacement list. Both bases use the same ordered candidate/orbit search.
Counts and finite-step errors can differ; the search does not guarantee
globally optimal conditioning. With `symmetry.saved_basis: true`, `symmetry.reduction: on` records
its basis and transformations in `fractional_basis.json` for either candidate
basis. By default this file is neither generated nor read. Existing cache files
are left untouched when the setting is disabled.

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
direction was selected in a fractional basis. Optional raw text exports preserve these
directional responses, not a manually selected subset of Cartesian slices.
With `symmetry.reduction: on`, symmetry may supply the negative sign; the text records this.

## Symmetry Denoising and Raman Tensors

For `atoms` and `on`, denoising fits the measured directional responses
within each representative's site-allowed, reciprocal tensor subspace before
transfer to equivalent atoms. For `off`, it projects the Cartesian field over
the complete space group, including equivalent-atom relationships. The latter
uses an equal-weight Cartesian norm, not an observation-space fit for
nonorthogonal fractional directions.

Reciprocity means $\partial\epsilon_{\alpha\beta}/\partial r_i
=\partial\epsilon_{\beta\alpha}/\partial r_i$, without complex conjugation;
do not impose it on general magneto-optical responses. No acoustic sum rule
is imposed, and symmetry-allowed numerical errors can remain.

Following the legacy derivative-file convention, optional stored atomic responses
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

The volume factor is applied once during in-memory mode contraction. It is also
included in optional atomic exports and must not be applied again when importing them. The
eigenvectors are mass-weighted and normalized. These are Placzek-response
derivatives, not a full resonant scattering calculation.

## Collective Mode Displacements

As an alternative to atomic displacements, use `CONTCAR`, Gamma-point
`band.yaml`, `irreps.yaml`, and the `symmetry` file to select Raman-active modes:

Add these settings to the same `input.yaml` file:

```yaml
modes:
  q: 0.01 # normal-coordinate amplitude, Angstrom sqrt(amu)
  indices: auto # select Raman-active modes automatically
  minimum_frequency: 1.0 # exclude modes at or below this cm^-1
  displacements_directory: mode_displacements
  results_directory: mode_results
```

These are the defaults. Set `modes.indices: [6, 7]`, for example, to select specific
1-based mode indices instead of applying the activity and frequency filters.
The directory settings control both generation and subsequent tensor extraction.

```bash
spectropy mode-displacements
# Run VASP externally in mode_displacements/ra_POSCAR_mode*/.
spectropy mode-tensors
```

Each mode, including each degenerate partner, receives its own explicit
positive and negative displacement. Automatic selection includes Raman-active
modes across polarizations with frequencies above `modes.minimum_frequency`
(default 1 cm⁻¹).

$$
\delta r_i(\mu)=\pm Q_k\frac{e_i^k(\mu)}{\sqrt{m_\mu}},
\qquad
R(k)\approx\frac{V_{\rm cell}}{4\pi}
\frac{\epsilon(+Q_k)-\epsilon(-Q_k)}{2Q_k}.
$$

`modes.q` is in Å$\sqrt{\mathrm{amu}}$, not a maximum atomic displacement;
it is independent of the atomic `displacements.amplitude` setting.
The generated directory contains `mode_manifest.json`, `mode_index`,
`POSCAR_mode*`, and `ra_POSCAR_mode*/POSCAR`. Supplied `INCAR`, `KPOINTS`, and
`POTCAR` files are symlinked automatically into every mode-displacement directory;
no separate staging script is generated or needed. Calculations are run externally.

`mode-tensors` validates displaced geometries and requires both signs, central
differences, and `polarization.mode: specific`. It writes
`raw_mode_tensors_<E>eV.npz`, `Raman_tensor_<E>eV`, and
`Raman_intensity_complex_<E>eV`. Each mode tensor is extracted directly from
its own pair; atomic denoising is a separate workflow.

## Publication-Quality Plotting

`spectropy plot` generates Raman spectra plots with Lorentzian or Gaussian
broadening and labels prominent peaks with their available mode symmetries.
It scans the working directory tree, using the selected unprefixed intensities.
PNG plots are normalized for display; numerical broadened spectra are not.
LaTeX is used when installed, with a Matplotlib-font fallback.

## Optional Command-line Overrides

The plain commands above are sufficient when settings are in `input.yaml`.
For one-off changes, an explicit CLI option overrides the corresponding setting
for that invocation only; it does not rewrite `input.yaml`.

**Precedence: explicit CLI option → `input.yaml` setting → built-in default.**

| Command | Optional overrides → `input.yaml` setting |
| --- | --- |
| `spectropy displacements` | `--symmetry` → `symmetry.reduction`; `--basis` → `displacements.basis`. |
| `spectropy tensors` | `--difference` → `displacements.difference`; `--equilibrium-file` → `displacements.equilibrium_file`; `--saved-basis` / `--no-saved-basis` → `symmetry.saved_basis`. |
| `spectropy mode-displacements` | `--q` → `modes.q`; `--modes` → `modes.indices`; `--minimum-frequency` → `modes.minimum_frequency`; `--output` → `modes.displacements_directory`. |
| `spectropy mode-tensors` | `--directory` → `modes.displacements_directory`; `--output` → `modes.results_directory`; `--energies` → `spectrum.energies`. |

All subcommands accept `--workdir PATH` (or `-C PATH`, default `.`) to choose
the directory containing `input.yaml`, and `--help` to list options. The top-level
command accepts `--version`. `--equilibrium-vasprun` is an alias for
`--equilibrium-file`. `spectropy symmetry` writes equivalent-atom mappings to
`symmetry_operation_matrices`; `spectrum` reads final mode tensors and `plot` reads intensities; both use `input.yaml`
without additional scientific CLI options.

## Examples

| Example | Starting point |
| --- | --- |
| [MoS2_vasp](examples/MoS2_vasp/README.md) | Minimal Cartesian 0.01 Å benchmark with optical inputs, phonons, three displaced structures, XML outputs, and a static Raman spectrum. |
| [MoS2_qe — WIP](examples/MoS2_qe/README.md) | Documentation starter only; the QE interface and a runnable dataset are not included. |

## Python Package

The implementation lives under `src/spectropy`. Import helpers from the package,
for example:

```python
from spectropy.config import read_settings
from spectropy.tensors import run_tensors

settings = read_settings("input.yaml")
run_tensors()
```

The `spectropy` command and `python -m spectropy` provide the same CLI.
Older flat imports such as `from spectropy_config import read_settings` must be
updated to `from spectropy.config import read_settings`; no duplicate flat-module
wrappers are installed. Calculation settings, output formats, and numerical
algorithms are unchanged by the package move.
