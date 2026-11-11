# Change notes

## Unreleased — supported workflow simplification

- Supplied `INCAR`, `KPOINTS`, and `POTCAR` inputs are automatically symlinked
  into atomic, equilibrium, and collective-mode calculation directories.
  Existing inputs are never replaced; missing inputs produce a warning.
- `fractional_basis.json` generation is opt-in with `saved_basis: true`.
  Default reconstruction does not need or read it. No existing cache is deleted.
- Collective-mode generation no longer writes the redundant `runRA` input-copy
  script; automatic symlinks replace that staging step.
- Settings for displacement selection/basis, mode generation, and optional
  reference/basis/output paths can be kept in `input`. Explicit CLI options
  override the file for one invocation; omitted options use the file's settings.
- Geometry uses named `incident_polarization`, `scattered_polarization`, and
  `surface_normal` entries; examples and newly written settings use this format.
- Atomic `full`, `atoms`, and `minimal` modes use one Cartesian/fractional
  direction selector and preserve raw directional responses.
- Atomic denoising remains optional; its output no longer includes residual reports.
- Collective-mode generation explicitly enumerates both signs of every selected
  Raman-active mode. Tensor extraction uses each directly computed pair.
- Removed mode partner/sign reduction, mode-level denoising, symmetry-check
  commands, residual exports/plots, symbolic forms, and polarization-angle scans.
- Mode manifests now use schema 2. Regenerate the plan for the current extractor;
  existing calculation outputs are not changed or deleted.
- Public line-width settings use `broadening_hwhm` in cm^-1, default 1.0.
- VASP optical input accepts nine-component or standard six-component
  `vasprun.xml`. Missing transpose entries are completed by reciprocity with a
  warning; printed components remain independent. XML validation remains strict,
  and there is no OUTCAR fallback.

These changes do not launch DFT calculations or establish convergence.
