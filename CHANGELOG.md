# Change notes

## Unreleased — supported workflow simplification

- Modules now live in the `spectropy` package under `src/spectropy`, with names
  such as `config`, `displacements`, `tensors`, and `spectrum`. CLI commands and
  file formats are unchanged; `python -m spectropy` is also supported. Python
  consumers must use package-qualified imports. Removed only unused wrappers,
  one unused import, and unreachable XML-reader fallback checks.
- Raw directional derivatives now use labeled text, `raw_directional_derivative_<E>`,
  with real/imaginary 3×3 matrices and direction/sign metadata, instead of JSON.
  The obsolete JSON-to-Cartesian reconstruction helper is removed; spectra
  continue to consume only final mode tensors.
- The VASP tutorial enables raw directional derivatives, final denoised atomic
  derivatives, and sampled dielectric matrices. The per-atom output switch also
  exports `dielectric_tensor_<E>`; the default remains mode tensors only.
- `spectrum.energies` accepts an inclusive `{start, stop, step}` range as well
  as an explicit list. Filename collisions at two decimal places are rejected.
- The verified VASP example is now `examples/MoS2_vasp`, with symmetry-on/off
  comparison instructions and a generated static spectrum. A separate
  `MoS2_qe` documentation starter is explicitly marked WIP.
- Settings are grouped under `polarization`, `spectrum`, `displacements`,
  `symmetry`, `broadening`, `modes`, and `output`. Flat keys are no longer read.
  The fixed XML source no longer needs a setting or a mode-tensor CLI flag.
- `spectropy tensors` replaces `spectropy derivatives`: reconstruct, optionally
  denoise, and contract atomic responses in memory. By default only the final
  `Raman_tensor_<E>eV` analysis files are written. `spectropy spectrum` reads
  those full-precision tensors to calculate intensities without atomic files.
- `output.per_atom.enabled` optionally saves the final `epsilon_derivative_<E>`
  field. Its nested `raw_directional_derivatives` option enables measured-basis
  text for debugging. No duplicate corrected field or raw/corrected spectrum
  copies are generated; existing outputs are not deleted.
- SpectroPy settings now use `input.yaml`: native YAML lists, booleans, and
  `#` comments replace the positional/text parser. All stages and examples use
  the new filename; explicit CLI settings still override file values.
  Duplicate/unknown keys and invalid values are rejected. Existing calculation
  inputs are not migrated automatically.
- Advanced usage is in `docs.md`; the README retains the basic workflow.
- Displacement selection is now `symmetry.reduction: on` (minimal), `off` (full), or
  `atoms`, with `--symmetry` as the CLI override. The old `displacement_mode`
  key and `--mode` option are removed. The selection algorithms and existing
  displacement-record format are unchanged.
- `ra_pos_equilibrium` is generated and staged only for
  `displacements.difference: forward`. Central generation leaves any existing
  equilibrium directory untouched; forward-reference documentation is in
  advanced usage rather than the basic workflow.
- Supplied `INCAR`, `KPOINTS`, and `POTCAR` inputs are automatically symlinked
  into atomic, equilibrium, and collective-mode calculation directories.
  Existing inputs are never replaced; missing inputs produce a warning.
- `fractional_basis.json` generation is opt-in with `symmetry.saved_basis: true`.
  Default reconstruction does not need or read it. No existing cache is deleted.
- Collective-mode generation no longer writes the redundant `runRA` input-copy
  script; automatic symlinks replace that staging step.
- Settings for displacement selection/basis, mode generation, and optional
  reference/basis/output paths can be kept in `input.yaml`. Explicit CLI options
  override the file for one invocation; omitted options use the file's settings.
- Geometry uses `incident`, `scattered`, and `surface_normal` entries in the
  `polarization` section; examples and newly written settings use this format.
- Atomic `off`, `atoms`, and `on` selections use one Cartesian/fractional
  direction selector; raw directional responses can be exported for debugging.
- Atomic denoising remains optional; its output no longer includes residual reports.
- Collective-mode generation explicitly enumerates both signs of every selected
  Raman-active mode. Tensor extraction uses each directly computed pair.
- Removed mode partner/sign reduction, mode-level denoising, symmetry-check
  commands, residual exports/plots, symbolic forms, and polarization-angle scans.
- Mode manifests now use schema 2. Regenerate the plan for the current extractor;
  existing calculation outputs are not changed or deleted.
- Public line-width settings use `broadening.hwhm` in cm^-1, default 1.0.
- VASP optical input accepts nine-component or standard six-component
  `vasprun.xml`. Missing transpose entries are completed by reciprocity with a
  warning; printed components remain independent. XML validation remains strict,
  and there is no OUTCAR fallback.

These changes do not launch DFT calculations or establish convergence.
