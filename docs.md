# SpectroPy Documentation

> **Work in progress:** This reference describes the current implementation and is still being refined. For a first calculation, start with the [README](README.md) and the [MoS₂ tutorial](examples/MoS2_vasp/README.md).

SpectroPy prepares finite displacements and converts dielectric responses into
Gamma-point Raman tensors and polarization-resolved spectra. VASP performs the
electronic calculations; Phonopy supplies crystal symmetry, phonon frequencies,
and eigenvectors. SpectroPy does not submit or manage VASP jobs.

## Contents

1. [Calculation model and workflow](#calculation-model-and-workflow)
2. [Installation and dependencies](#installation-and-dependencies)
3. [Required scientific inputs](#required-scientific-inputs)
4. [Configuration reference](#configuration-reference)
5. [Atomic displacement generation](#atomic-displacement-generation)
6. [Importing Phonopy displacements](#importing-phonopy-displacements)
7. [Dielectric data and finite differences](#dielectric-data-and-finite-differences)
8. [Reconstruction, denoising, and mode contraction](#reconstruction-denoising-and-mode-contraction)
9. [Spectra and plotting](#spectra-and-plotting)
10. [Collective mode workflow](#collective-mode-workflow)
11. [Files and output formats](#files-and-output-formats)
12. [Command-line reference](#command-line-reference)
13. [Mode visualization](#mode-visualization)
14. [Validation and troubleshooting](#validation-and-troubleshooting)
15. [Python package and examples](#python-package-and-examples)

## Calculation model and workflow

Two routes produce the same mode-tensor file format:

```text
Atomic route
  reference structure + symmetry
    → spectropy displacements → VASP (external)
    → spectropy tensors ← Gamma phonons
        directional responses → Cartesian atomic field
        → optional denoising/reconstruction → mode tensors

Collective mode route
  reference structure + Gamma phonons + symmetry
    → spectropy mode-displacements → VASP (external)
    → spectropy mode-tensors → mode tensors

Either route
  mode tensors → spectropy spectrum → spectropy plot
```

The atomic route reconstructs a response field reusable for all modes of the
reference cell. The collective route displaces atoms together along selected
normal modes. These are alternatives, not consecutive stages.

| Quantity | Convention |
| --- | --- |
| Atomic positions and displacement lengths | Å |
| Stored displacement vectors | Fractional finite vectors in `displacements.dat` |
| Dielectric tensor indices | Cartesian; dimensionless complex response |
| Final atomic-response indices | Cartesian |
| Laser energy | eV |
| Phonon frequencies | THz and/or cm⁻¹, as labeled |
| Normal-coordinate amplitude | Å√amu |
| Broadening | HWHM in cm⁻¹ |

## Installation and dependencies

Use Python 3.9 or newer. From a checkout:

```bash
pip install -e '.[all]'
spectropy --help
```

| Dependency | Role |
| --- | --- |
| NumPy, PyYAML | Core numerical operations and settings; installed automatically |
| Matplotlib (`.[plot]`) | Spectrum plotting |
| Phonopy (`.[phonopy]`) | Phonon preparation and displacement import |
| `.[all]` | Installs plotting, Phonopy, and spglib dependencies |
| VASP | External dielectric/force calculations; separately installed and licensed |
| VESTA or VMD | Optional viewing of phonon animations/arrows |

After updating a checkout, reinstall if entry points or dependencies changed.
`python -m spectropy` provides the same commands as `spectropy`.

## Required scientific inputs

### Inputs by stage

| Stage | Inputs |
| --- | --- |
| Atomic displacement preparation | `CONTCAR`, `symmetry`, `input.yaml`; supplied `INCAR`, `KPOINTS`, `POTCAR` for staging |
| Phonopy import | The same inputs plus the selected `phonopy_disp.yaml` |
| Atomic tensor calculation | Reference structure, displacement records, symmetry, complete optical XMLs, `band.yaml`; `irreps.yaml` for labels |
| Collective mode preparation | Reference structure, `band.yaml`, `irreps.yaml`, symmetry, settings and VASP inputs |
| Spectrum calculation | `Raman_tensor_<E>eV`, settings with polarization vectors |
| Plotting | Intensity files and broadening settings |

### Preparing Phonopy results

From a completed force calculation, use the reference structure and the same
supercell dimensions that were used to generate its displacements:

```bash
phonon_dim="4 4 1" # Replace with the dimensions used for your force calculations.
phonopy -f hf_POSCAR-*/vasprun.xml
phonopy -c CONTCAR --symmetry > symmetry
phonopy -c CONTCAR --dim="$phonon_dim" \
  --band="0 0 0 0 0 0" --band-points=1 --eigvecs
phonopy -c CONTCAR --dim="$phonon_dim" --irreps="0 0 0"
```

Run force collection with the original `phonopy_disp.yaml` present. The shell
glob assumes zero-padded directory names whose sorted order matches the YAML;
otherwise list the XML paths explicitly in displacement-record order.

Phonopy reads `FORCE_SETS` and builds the force constants internally for the
band and irrep calculations. Our workflow does not require a saved
`FORCE_CONSTANTS` file or `--readfc`. These XMLs must contain converged forces
for the original Phonopy displacements, not an unrelated Raman displacement set.
See the [Phonopy force-calculation instructions](https://phonopy.github.io/phonopy/vasp.html).

To subtract equilibrium forces in a **combined force/Raman workflow**, use
`phonopy --fz` instead of `phonopy -f`, listing the undisplaced equilibrium XML **first**, followed by
every displaced XML in `phonopy_disp.yaml` order. This subtracts equilibrium
forces before writing `FORCE_SETS`; the band and irrep commands remain the same.

`-c CONTCAR` specifies the **reference cell**, not the enlarged force-calculation
`SPOSCAR`. `--dim` specifies the enlargement used for those force constants;
it is **not** the Raman displacement-cell size. For example, `4 4 1` uses
48-atom force calculations but produces nine Gamma modes for the three-atom
MoS₂ reference cell. Use `1 1 1` only if the force calculation used no enlargement.

The Gamma-to-Gamma path writes a Gamma-only `band.yaml`;
SpectroPy reads the first q-point. `--eigvecs` is essential; band labels and
a plotted dispersion are not. `--qpoints="0 0 0"` writes `qpoints.yaml`, which
is **not** a substitute for `band.yaml`. See Phonopy's
[output formats](https://phonopy.github.io/phonopy/output-files.html) and
[band settings](https://phonopy.github.io/phonopy/setting-tags.html#band-structure-tags).

Copy `band.yaml`, `symmetry`, and (if generated) `irreps.yaml` alongside the
same `CONTCAR` in the Raman directory. SpectroPy does **not** need `FORCE_SETS`,
`FORCE_CONSTANTS`, or the force-job XMLs for the atomic workflow shown here.
`phonopy_disp.yaml` is only needed when explicitly importing Phonopy's
displacements, not when SpectroPy generates its own displacements.

### Registering phonons against the response structure

| Registration check | Requirement |
| --- | --- |
| Atom correspondence | Species and periodic positions; unique common origin shift allowed |
| Lattice | Same Cartesian frame; component differences ≤10⁻⁶ Å |
| Positions | Matching tolerance 10⁻³ Å |
| Alignment | Eigenvectors and masses reordered in memory; source files unchanged |
| Labels | Read from `irreps.yaml`; optional for atomic contraction, needed for automatic mode selection |

`band.yaml` supplies eigenvectors, masses, and frequencies. Ambiguous atom matches fail.

## Configuration reference

### Working configuration

This atomic-route example uses minimal Cartesian displacements and retains
optional atomic/debug outputs. It is a coherent example, not the default for
every setting.

```yaml
polarization:
  incident: [1, 0, 0]
  scattered: [1, 0, 0]
  surface_normal: z
  mode: specific
displacements:
  basis: cartesian
  amplitude: 0.01
  difference: central
  equilibrium: true
symmetry:
  reduction: on
  denoise: true
spectrum:
  energies: {start: 0.0, stop: 3.0, step: 0.1}
  cutoff: 0
  temperature: 298
broadening:
  type: lorentzian
  hwhm: 1
output:
  per_atom:
    enabled: true
    raw_directional_derivatives: true
```

Mode-specific settings are documented under [Collective mode workflow](#collective-mode-workflow);
export settings are documented under [Files and output formats](#files-and-output-formats).

Settings are grouped YAML mappings in `input.yaml`; their order does not matter.
The table uses dotted paths (for example, `symmetry.denoise` means `denoise`
inside `symmetry`), not literal dotted YAML keys. See the basic example.
Use lists for vectors and multiple values, `true`/`false` for booleans, and
`#` for comments. Duplicate and unknown keys are rejected.
`symmetry.reduction: on` and `symmetry.reduction: off` are string choices, not booleans; quotes are
optional in SpectroPy. The loader uses only `true` and `false` as boolean words.

### Polarization

| Setting | Default | Meaning |
| --- | --- | --- |
| `polarization.incident` | Required for spectra | Three Cartesian components of the incoming polarization vector. |
| `polarization.scattered` | Required for spectra | Three Cartesian components of the detected polarization vector. |
| `polarization.surface_normal` | Unset | Optional geometry label: `x`, `y`, or `z`. Does not alter the explicitly specified polarization vectors or intensities. |
| `polarization.mode` | `specific` | Only `specific` is supported: calculate the intensity for the supplied incident and scattered vectors. |

### Laser energies and spectral factors

| Setting | Default | Meaning |
| --- | --- | --- |
| `spectrum.energies` | `0.00` | A scalar, list, or inclusive `{start, stop, step}` energy range in eV; zero is the static limit. |
| `spectrum.cutoff` | Unset | Apply the inverse phonon-frequency factor only to positive modes at or above this threshold (cm⁻¹). Unset disables the denominator; modes below the cutoff retain a denominator factor of 1, not zero intensity. |
| `spectrum.temperature` | `298` | Temperature in kelvin for the independent Stokes Bose factor, applied to all positive-frequency modes regardless of cutoff. At 0 K this factor is 1. Zero/negative-frequency modes receive neither correction. |

### Atomic displacements

| Setting | Default | Meaning |
| --- | --- | --- |
| `displacements.basis` | `cartesian` | SpectroPy generation: `cartesian` or `fractional` candidate directions. Omit for Phonopy imports. |
| `displacements.amplitude` | `0.03` | SpectroPy generation: physical displacement length in Å, in either basis. Omit for Phonopy imports. |
| `displacements.phonopy_file` | Unset | Import Phonopy's recorded vectors instead of generating a displacement set. |
| `displacements.difference` | `central` | SpectroPy **generation only**: `central` prepares opposite signs (reducing redundant signs with symmetry); `forward` prepares one sign and equilibrium. Omit for Phonopy imports. Tensor evaluation infers the formula from the available records. |
| `displacements.equilibrium_file` | Automatic | Optional equilibrium XML path for forward differences; otherwise `displacements/ra_pos_equilibrium/vasprun.xml`. |
| `displacements.equilibrium` | `true` | Generate `ra_pos_equilibrium` alongside the displaced structures, including central runs. Set `false` to skip the additional directory; existing data are not deleted. |

### Symmetry

| Setting | Default | Meaning |
| --- | --- | --- |
| `symmetry.reduction` | `off` | `on` (minimal), `off` (full), or `atoms` (atoms only); see the selection table below. |
| `symmetry.denoise` | `false` | `true` enables atomic symmetry and reciprocity denoising. |
| `symmetry.saved_basis` | `false` | For `on`, generate `fractional_basis.json` and reuse it after validation; otherwise reconstruct without a cache. |

### Line broadening

| Setting | Default | Meaning |
| --- | --- | --- |
| `broadening.hwhm` | `1.0` | Half-width at half-maximum in cm⁻¹; FWHM is twice this value. |
| `broadening.type` | `lorentzian` | `lorentzian` or `gaussian` (also `l` or `g`). |

| Configuration rule | Behavior |
| --- | --- |
| Amplitude list | `off`/`atoms` accept three axis lengths; `on` requires one common length |
| Paths | Relative to the working directory |
| Existing calculations | Editing YAML does not alter structures or results |
| Changed displacement settings | Generate a fresh dataset; analysis uses recorded vectors |
| Older input formats | Positional and flat settings must be rewritten using the grouped schema |
| Dielectric source | `vasprun.xml` only |

## Atomic displacement generation

| `symmetry.reduction` | Atoms | Directions |
| --- | --- | --- |
| `off` | Every atom | Three selected-basis axes. |
| `atoms` | One representative per atomic orbit | Three selected-basis axes. |
| `on` | One representative per atomic orbit | Site-symmetry-reduced axis or mixed directions. |

| Selection detail | Behavior |
| --- | --- |
| Reduction vs denoising | Independent: `off` can still use `denoise: true` |
| Cartesian vs fractional | Same candidate/orbit search; direction counts and finite-step errors can differ |
| Optimality | Globally optimal conditioning is not guaranteed |
| Default source | SpectroPy selects directions; Phonopy supplies symmetry |
| Optional cache | `saved_basis: true` with reduction `on` writes/reads `fractional_basis.json`; disabled caches are ignored |

### Staged calculation layout

```text
calculation/
├── CONTCAR
├── input.yaml
├── symmetry
├── INCAR, KPOINTS, POTCAR
├── displacements.dat
└── displacements/
    ├── ra_pos_atom…/
    │   ├── POSCAR
    │   └── INCAR, KPOINTS, POTCAR → supplied inputs
    └── ra_pos_equilibrium/
        ├── POSCAR
        └── INCAR, KPOINTS, POTCAR → supplied inputs
```

Input symlinks are created for supplied files. Missing inputs produce a warning;
matching links are retained, and conflicting files or links are never overwritten.
Supply any missing inputs before running VASP.

### Equilibrium calculations

| Equilibrium setting | Behavior |
| --- | --- |
| Default | Generate `displacements/ra_pos_equilibrium` for central and forward sets |
| `displacements.equilibrium: false` | Skip extra generation; existing data remain |
| `displacements.equilibrium_file` | Override the equilibrium XML path |
| Required for analysis | Only when neither a measured nor symmetry-derived opposite sign is available |

Use matching optical settings for equilibrium and displaced calculations.

## Importing Phonopy displacements

To use an existing VASP `phonopy_disp.yaml` instead of selecting new directions:

```yaml
displacements:
  phonopy_file: phonopy_disp.yaml
symmetry:
  reduction: on
  denoise: true
```

Run `spectropy displacements` in a fresh directory with matching structure,
symmetry and VASP inputs. Requires `spectropy[phonopy]`. After VASP, use the
ordinary `spectropy tensors` pipeline; the selector is bypassed, not the
reconstruction or denoising.


| Item | Import behavior |
| --- | --- |
| Reference | `CONTCAR` must match the YAML **supercell**, including atom order and Cartesian frame—not its smaller unit cell. |
| Directions and lengths | Taken from the YAML unchanged. Omit `displacements.basis`, `amplitude`, and `difference`; specifying them with `phonopy_file` is an error. |
| Metadata | Writes the usual `displacements.dat` with fractional finite vectors; reconstruction still produces Cartesian derivatives. |
| Calculation folders | Uses the usual `displacements/ra_pos_atom*` layout and links supplied `INCAR`, `KPOINTS`, and `POTCAR`. |
| Record mapping | `phonopy_displacement_map.tsv` maps each YAML displacement record (1-based, in file order) to its SpectroPy directory and Cartesian vector. |
| Existing data | Refuses to overwrite a displacement set. It does **not** automatically copy/link old VASP results. |

| Dataset check | Requirement |
| --- | --- |
| Supported records | Single-atom VASP/Å displacements; zero vectors and ± partners recognized |
| Rejected records | Multi-atom snapshots, other calculator units, repeated amplitudes along one direction, unequal opposite steps |
| Equilibrium | Added by default; does not override central-pair precedence |
| Symmetry | Must describe the same electronic/magnetic state |
| Reusing results | Match atom, signed vector, amplitude, geometry and electronic settings—not directory numbers |
| Optical response | Force-only XMLs are insufficient |

Import preserves YAML order and adds no nonzero displacements.

## Dielectric data and finite differences

### XML compatibility

| XML contents | Reader behavior |
| --- | --- |
| Nine dielectric components | Preserve each independently |
| Six components | Warn and fill transpose entries by reciprocity, without complex conjugation |
| Malformed XML, missing diagonal or entire off-diagonal pair | Fail |
| OUTCAR | Not used |

Inferred transpose entries are not measurements and cannot recover a nonreciprocal response.

### Directional derivatives

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

### Formula inference and sign precedence

For either displacement source, `spectropy tensors` classifies the actual vectors
per atom and direction (including mixed directions, not just Cartesian axes):

| Available observations | Inferred derivative |
| --- | --- |
| Equal-length measured $+h\mathbf u$ and $-h\mathbf u$ | Central: $[\epsilon(+h\mathbf u)-\epsilon(-h\mathbf u)]/(2h)$. |
| One signed displacement and a valid site sign-reversal operation | Central, with the missing response supplied by symmetry, even when equilibrium is present. |
| One signed displacement and equilibrium, without a valid sign reversal | Forward along the recorded signed direction: $[\epsilon(h\mathbf u)-\epsilon(0)]/h$. |
| Neither a partner nor an equilibrium response nor a valid sign reversal | Error; no missing response is assumed zero. |

Precedence: **measured pair → symmetry partner → equilibrium → error**.
Each direction is classified separately; raw exports record the choice.
Negative-only records use their signed direction. Zero vectors are equilibrium,
never derivatives.

## Reconstruction, denoising, and mode contraction

| Stage | Operation |
| --- | --- |
| Cartesian reconstruction | Solve directional observations for the rank-three atomic field |
| Denoising: `on`/`atoms` | Fit measured responses in each representative's site-allowed reciprocal subspace, then transfer |
| Denoising: `off` | Project the full Cartesian field over the space group, including equivalent atoms |
| Atom transfer | Rotate all three tensor indices, including the displacement index |

Full-field projection uses an equal-weight Cartesian norm, not an
observation-space fit for nonorthogonal directions. Reciprocity means
$D_{i,\alpha\beta}=D_{i,\beta\alpha}$, **without conjugation**; it is not
appropriate for general magneto-optical responses. No acoustic sum rule is
imposed. Symmetry-allowed errors can remain.

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

## Spectra and plotting

Polarization selects the contraction of the complex mode tensor:

$$
I_k^{\mathrm{reduced}}=|\mathbf g_s^{\mathsf T}R(k)\mathbf g_i|^2.
$$

Only the specified polarization geometry is implemented. The surface-normal
label does not rotate or replace the supplied polarization vectors.

For positive-frequency modes, the Stokes population factor is
$ n_k(T)+1=1/[1-\exp(-\hbar\omega_k/k_B T)] $; at zero temperature it is one.
The independent `spectrum.cutoff` setting controls the frequency denominator,
not this population factor and not a hard intensity filter. Refer to the
configuration table for behavior below the cutoff. Intensities are relative,
not calibrated absolute scattering cross sections.

### Energy selection

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

### Broadening and figures

`spectropy plot` generates Raman spectra plots with Lorentzian or Gaussian
broadening and labels prominent peaks with their available mode symmetries.
It scans the working directory tree, using the selected unprefixed intensities.
PNG plots are normalized for display; numerical broadened spectra are not.
LaTeX is used when installed, with a Matplotlib-font fallback.

### What to rerun after changing settings

| Change | Repeat |
| --- | --- |
| Displacement size, basis, selection, or reference structure | Prepare a fresh dataset and run its electronic calculations |
| Laser energies or denoising | `tensors`, `spectrum`, `plot` |
| Polarization, temperature, or frequency cutoff | `spectrum`, `plot` |
| Broadening width or line shape | `plot` |
| Phonon eigenvectors/masses | Tensor contraction and downstream stages |


## Collective mode workflow

> **Not yet tested:** The collective mode workflow has not yet been validated end to end. Treat its results as unverified until independently checked against a trusted calculation.

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

| Mode setting | Meaning |
| --- | --- |
| `modes.q` | Normal-coordinate amplitude in Å√amu, not maximum atomic motion |
| `modes.indices: auto` | All Raman-active partners above `minimum_frequency`, across polarizations |
| `modes.indices: [6, 7]` | Explicit 1-based indices; bypass automatic activity/frequency filtering |
| Directory settings | Used by both generation and extraction |

```bash
spectropy mode-displacements
# Run VASP externally in mode_displacements/ra_POSCAR_mode*/.
spectropy mode-tensors
# Use the same settings for the ordinary spectrum stage (default output folder).
ln -s ../input.yaml mode_results/input.yaml
spectropy spectrum --workdir mode_results
```

Each selected mode—including every degenerate partner—gets explicit ± displacements.

$$
\delta r_i(\mu)=\pm Q_k\frac{e_i^k(\mu)}{\sqrt{m_\mu}},
\qquad
R(k)\approx\frac{V_{\rm cell}}{4\pi}
\frac{\epsilon(+Q_k)-\epsilon(-Q_k)}{2Q_k}.
$$

| Mode-workflow output | Contents |
| --- | --- |
| `mode_manifest.json`, `mode_index` | Mode/displacement metadata |
| `POSCAR_mode*`, `ra_POSCAR_mode*/POSCAR` | Collective displacements |
| Calculation directories | Symlinks to supplied INCAR, KPOINTS and POTCAR |
| `raw_mode_tensors_<E>eV.npz` | Direct mode-tensor data |
| `Raman_tensor_<E>eV` | Spectrum-compatible tensors with original indices and available labels |

`mode-tensors` validates geometries and requires both signs. It does not require
polarization settings or apply atomic denoising. Run `spectrum` in the results
directory with the desired polarization and spectral settings.

## Files and output formats

### Displacement metadata and auxiliary files

| File | Purpose |
| --- | --- |
| `fractional_basis.json` | Optional JSON reconstruction cache, generated and used only when both `symmetry.saved_basis: true` and `symmetry.reduction: on` are set. |
| `atomic_displacement` | Fortran-format export of the displacement list; generated, but not read by SpectroPy. |
| `displacements/pos_atom*` | Standalone POSCAR exports duplicating the calculation-directory structures; not read by SpectroPy. |

Only `displacements.dat` and the per-calculation `POSCAR` files are essential
generated inputs. Keep the original `CONTCAR` and symmetry/phonon inputs too.
The tensor stage registers phonon atoms against the reconstruction's `CONTCAR`.
`irreps.yaml` supplies labels and can be omitted for atomic-tensor calculations.

### Tensor and dielectric exports

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

| Export detail | Convention |
| --- | --- |
| Final atomic field | Denoised only when `symmetry.denoise: true`; no duplicate corrected file |
| Raw directional text | Measured directions, steps, sign provenance, and real/imaginary 3×3 responses—not Cartesian slices |
| Sampled dielectric tensors | Dimensionless, Cartesian, nearest-grid XML values before denoising; signed vectors in Å |
| Equilibrium sample | Included when used for forward differentiation |
| Six-component sources | Usual transpose-completion warning still applies |
| Numerical precision | Full-precision tensor text; no extra XML reads for these exports |
| Spectrum input | Final mode tensors only; atomic/debug exports are optional |
| Old exports | Not automatically deleted or selected |

### Mode tensor and intensity records

| File | Record structure |
| --- | --- |
| `Raman_tensor_<E>eV` | Per-mode header: index, frequency in THz, frequency in cm⁻¹, irrep; followed by three rows of complex `real±imagj` values |
| `epsilon_derivative_<E>` | Separate real/imaginary sections; each atom has a 3×9 block containing the three Cartesian 3×3 derivative slices side by side |
| `Raman_intensity_complex_<E>eV` | Mode frequency in cm⁻¹, real nonnegative intensity, and available irrep label; “complex” refers to the input tensor, not a complex intensity |

Energy-qualified filenames use two decimal places. Do not round tensor data
for subsequent analysis; round only presentation tables.

## Command-line reference

**Precedence: explicit CLI option → YAML → default.** Overrides affect only the
current invocation; they do not rewrite `input.yaml`.


| Command | Optional overrides → `input.yaml` setting |
| --- | --- |
| `spectropy displacements` | `--symmetry` → `symmetry.reduction`; `--basis` → `displacements.basis`; `--phonopy` → `displacements.phonopy_file`. |
| `spectropy tensors` | `--equilibrium-file` → `displacements.equilibrium_file`; `--saved-basis` / `--no-saved-basis` → `symmetry.saved_basis`. The finite-difference formula is inferred, not overridden. |
| `spectropy mode-displacements` | `--q` → `modes.q`; `--modes` → `modes.indices`; `--minimum-frequency` → `modes.minimum_frequency`; `--output` → `modes.displacements_directory`. |
| `spectropy mode-tensors` | `--directory` → `modes.displacements_directory`; `--output` → `modes.results_directory`; `--energies` → `spectrum.energies`. |

All subcommands accept `--workdir PATH` (or `-C PATH`, default `.`) to choose
the directory containing `input.yaml`, and `--help` to list options. The top-level
command accepts `--version`. `--equilibrium-vasprun` is an alias for
`--equilibrium-file`. `spectropy symmetry` writes equivalent-atom mappings to
`symmetry_operation_matrices`; `spectrum` reads final mode tensors and `plot` reads intensities; both use `input.yaml`
without additional scientific CLI options.

## Mode visualization

This script generates files to visualize the atomic motion of each phonon mode in either **VMD** or **VESTA**.

For VESTA, open `CONTCAR` and save a project file (e.g. `template.vesta`) in
your working directory. From that directory, run the installed helper:

```bash
python -m spectropy.visualization
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

## Validation and troubleshooting

| Symptom | Check |
| --- | --- |
| Unknown/duplicate YAML key | Use the grouped schema and spellings in the configuration reference |
| Missing VASP input link | Supply the input before generation; inspect any conflicting file rather than overwriting it |
| Missing opposite displacement | Check recorded vectors and site symmetry; supply equilibrium only when a central partner cannot be obtained |
| Incomplete reconstruction coverage | Ensure representative directions and their symmetry images span the required displacement space |
| Phonon registration failure | Compare atom counts, species, periodic positions, and Cartesian lattice frames |
| XML parsing failure | Confirm complete, converged optical output; there is no OUTCAR fallback |
| Six-component XML warning | Transpose entries were inferred, not independently measured |
| Energy outside the XML grid | Request energies within the computed grid |
| Missing irrep labels | Provide a matching `irreps.yaml`; do not invent labels from frequency alone |
| Unexpectedly different plots | Check polarization, temperature, cutoff, broadening, normalization, and eigenvector basis |
| Large change after denoising | Check geometry, magnetic symmetry, numerical convergence and displacement dependence; projection alone does not establish accuracy |

Changing settings does not validate an old calculation under the new settings.
Keep structures, displacement metadata, XMLs and phonons together. For comparisons,
use a common eigenvector basis or account explicitly for mode signs and degenerate
partner rotations. A sign change of a mode tensor alone need not change intensity.

## Python package and examples

The implementation lives under `src/spectropy`. Import helpers from the package,
for example:

```python
from spectropy.config import read_settings
from spectropy.tensors import run_tensors

settings = read_settings("input.yaml")
run_tensors(settings=settings)  # Reads scientific files from the working directory.
```

Use package imports under `spectropy`; older flat-module imports are unsupported.

| Example | Starting point |
| --- | --- |
| [MoS2_vasp](examples/MoS2_vasp/README.md) | Minimal Cartesian 0.01 Å benchmark with optical inputs, phonons, three displaced structures, XML outputs, and a static Raman spectrum. |
| [MoS2_qe — WIP](examples/MoS2_qe/README.md) | Documentation starter only; the QE interface and a runnable dataset are not included. |
