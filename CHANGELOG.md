# Change notes

## Unreleased — supported workflow simplification

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
- VASP optical input is strictly nine-component `vasprun.xml`; there is no
  OUTCAR fallback or automatic transpose completion.

These changes do not launch DFT calculations or establish convergence.
