# SpectroPy Advanced Usage

See the [basic workflow](README.md#workflow--usage-guide) for a short introduction.

## Input Settings

Settings are a YAML mapping in `input.yaml`; their order does not matter.
Use lists for vectors and multiple values, `true`/`false` for booleans, and
`#` for comments. Duplicate and unknown keys are rejected.
`symmetry: on` and `symmetry: off` are string choices, not booleans; quotes are
optional in SpectroPy. The loader uses only `true` and `false` as boolean words.

| Setting | Default | Meaning |
| --- | --- | --- |
| `incident_polarization` | Required for spectra | Three Cartesian components of the incoming polarization vector. |
| `scattered_polarization` | Required for spectra | Three Cartesian components of the detected polarization vector. |
| `surface_normal` | Required for spectra | Surface-normal axis: `x`, `y`, or `z`. |
| `laser_energies` | `0.00` | One or more photon energies in eV; zero is the static limit. |
| `symmetry` | `off` | `on` (minimal), `off` (full), or `atoms` (atoms only); see the selection table below. |
| `displacement_basis` | `cartesian` | `cartesian` or `fractional` candidate directions. |
| `displacement_amplitude` | `0.03` | Physical displacement length in Å, in either basis. |
| `difference_scheme` | `central` | `central` or `forward` derivative analysis. Only `forward` generates an equilibrium calculation directory; the signed displacement set is unchanged. |
| `dielectric_source` | `vasprun.xml` | Only supported dielectric source. |
| `site_symmetry_projection` | `false` | `true` enables atomic symmetry and reciprocity denoising. |
| `polarization` | `specific` | `specific` or `average`; the latter sums polarization channels, not a general powder orientational integral. |
| `cutoff` | Unset | When set, apply the 298 K Bose/phonon-frequency factor to modes at or above this threshold (cm⁻¹). Otherwise use reduced intensities. |
| `broadening_hwhm` | `1.0` | Half-width at half-maximum in cm⁻¹; FWHM is twice this value. |
| `broadening_type` | `lorentzian` | `lorentzian` or `gaussian` (also `l` or `g`). |
| `equilibrium_file` | Automatic | Optional equilibrium XML path for forward differences; otherwise `displacements/ra_pos_equilibrium/vasprun.xml`. |
| `saved_basis` | `false` | For `on`, generate `fractional_basis.json` and reuse it after validation; otherwise reconstruct without a cache. |

`off` and `atoms` accept a list of three displacement lengths, one per selected
basis axis, e.g. `displacement_amplitude: [0.01, 0.02, 0.03]`; `on` requires
a common length. Accepted aliases are `laser_energy` and `frequency_cutoff`.
Legacy `broadening_fwhm` is converted to HWHM; do not specify both width keys.

The former text file named `input` is no longer read automatically. To migrate,
create `input.yaml` using named entries, replace space-separated vectors with
YAML lists, and replace `!` comments with `#`. Renaming alone is insufficient.
Replace `displacement_mode: minimal/full/atoms` with `symmetry: on/off/atoms`,
respectively. The CLI override is now `--symmetry`, not `--mode`.
If only the old file is present, SpectroPy reports the required conversion
instead of silently using defaults. Existing calculation files are not modified.

Paths are relative to the working directory. Changing `input.yaml` does not modify
existing structures or results: regenerate displacements in a fresh directory
to change their length, selection, or basis. Derivative analysis uses the
recorded displacement metadata, not a newly edited selection setting.

## Auxiliary Files and Optional Behavior

| File | Purpose |
| --- | --- |
| `fractional_basis.json` | Optional JSON reconstruction cache, generated and used only when both `saved_basis: true` and `symmetry: on` are set. |
| `atomic_displacement` | Fortran-format export of the displacement list; generated, but not read by SpectroPy. |
| `displacements/pos_atom*` | Standalone POSCAR exports duplicating the calculation-directory structures; not read by SpectroPy. |

Only `displacements.dat` and the per-calculation `POSCAR` files are essential
generated inputs. Keep the original `CONTCAR` and symmetry/phonon inputs too.
Atom registration prefers `ref_poscar.vasp` but can fall back to `CONTCAR`.
`irreps.yaml` supplies labels and can be omitted for atomic-tensor calculations.

For forward differences, set `difference_scheme: forward` before generating
displacements. This also creates `displacements/ra_pos_equilibrium/POSCAR`,
copied from the reference structure, with the same input symlinks. Run VASP
there with matching optical settings; analysis automatically reads its
`vasprun.xml`, or an existing response specified by `equilibrium_file`.
Central differences do not create or update this directory. Existing
equilibrium calculations are left untouched, not deleted.

Input symlinks are created for supplied files. Missing inputs produce a warning;
matching links are retained, and conflicting files or links are never overwritten.
Supply any missing inputs before running VASP.

Disabling `site_symmetry_projection` selects the raw result for the unprefixed
outputs. Existing corrected files are not deleted; only the active selection
changes when the analysis is rerun.

## Automating Energy Scans

Request multiple energies in `input.yaml`, then run the analysis once:

```yaml
laser_energies: [0.00, 0.10, 0.20, 0.30, 1.50, 1.96, 2.33]
```

```bash
spectropy derivatives
spectropy spectrum
spectropy plot
```

The nearest tabulated energy is selected without interpolation. Filenames use
two decimal places; choose energies that remain distinct at that precision.

## Displacement Selection and Derivatives

| `symmetry` | Atoms | Directions |
| --- | --- | --- |
| `off` | Every atom | Three selected-basis axes. |
| `atoms` | One representative per atomic orbit | Three selected-basis axes. |
| `on` | One representative per atomic orbit | Site-symmetry-reduced axis or mixed directions. |

`symmetry` controls displacement selection, independently of denoising:
`symmetry: off` still allows `site_symmetry_projection: true` when comparing
the full set against a symmetry-denoised result. The existing algorithm names
in `displacements.dat` metadata are retained so completed datasets remain readable.

All three support Cartesian or fractional candidates. SpectroPy generates
displacements itself; Phonopy supplies symmetry information, not the
displacement list. Both bases use the same ordered candidate/orbit search.
Counts and finite-step errors can differ; the search does not guarantee
globally optimal conditioning. With `saved_basis: true`, `symmetry: on` records
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
direction was selected in a fractional basis. Raw JSON records preserve these
directional responses, not a manually selected subset of Cartesian slices.
With `symmetry: on`, symmetry may supply the negative sign; the JSON records this.

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

## Collective Mode Displacements

As an alternative to atomic displacements, use `CONTCAR`, Gamma-point
`band.yaml`, `irreps.yaml`, and `symmetry` to select Raman-active modes:

Add these settings to the same `input.yaml` file:

```yaml
mode_q: 0.01 # normal-coordinate amplitude, Angstrom sqrt(amu)
mode_indices: auto # select Raman-active modes automatically
mode_minimum_frequency: 1.0 # exclude modes at or below this cm^-1
mode_displacements_directory: mode_displacements
mode_results_directory: mode_results
```

These are the defaults. Set `mode_indices: [6, 7]`, for example, to select specific
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
`POSCAR_mode*`, and `ra_POSCAR_mode*/POSCAR`. Supplied `INCAR`, `KPOINTS`, and
`POTCAR` files are symlinked automatically into every mode-displacement directory;
no separate staging script is generated or needed. Calculations are run externally.

`mode-tensors` validates displaced geometries and requires both signs, central
differences, and `polarization: specific`. It writes
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
| `spectropy displacements` | `--symmetry` → `symmetry`; `--basis` → `displacement_basis`. |
| `spectropy derivatives` | `--difference` → `difference_scheme`; `--equilibrium-file` → `equilibrium_file`; `--saved-basis` / `--no-saved-basis` → `saved_basis`. |
| `spectropy mode-displacements` | `--q` → `mode_q`; `--modes` → `mode_indices`; `--minimum-frequency` → `mode_minimum_frequency`; `--output` → `mode_displacements_directory`. |
| `spectropy mode-tensors` | `--directory` → `mode_displacements_directory`; `--output` → `mode_results_directory`; `--source` → `dielectric_source`; `--energies` → `laser_energies`. |

All subcommands accept `--workdir PATH` (or `-C PATH`, default `.`) to choose
the directory containing `input.yaml`, and `--help` to list options. The top-level
command accepts `--version`. `--equilibrium-vasprun` is an alias for
`--equilibrium-file`. `spectropy symmetry` writes equivalent-atom mappings to
`symmetry_operation_matrices`; `spectrum` and `plot` use their `input.yaml` settings
without additional scientific CLI options.

## Examples

| Example | Starting point |
| --- | --- |
| [MoS2_Tutorial](examples/MoS2_Tutorial/README.md) | Bundled derivative data: calculate tensors and spectra. |
| [MoS2_From_Displacements](examples/MoS2_From_Displacements/README.md) | Detailed earlier-stage walkthrough; structure, phonon, and raw dielectric bundles are placeholders, not yet supplied. |
| [MoS2_Mode_Displacements](examples/MoS2_Mode_Displacements/README.md) | Displace collective Raman-active modes instead of individual atoms. |
