# Change notes

## Unreleased — VASP workflow updates

- Displacement selection now uses independent `--mode full|atoms|minimal` and
  `--basis cartesian|fractional` options with the same direction-search algorithm.
  The old `--mode fractional` spelling has been removed.
- Raw atomic responses are saved as measured-direction JSON, not as a raw
  Cartesian tensor field. Denoising fits site symmetry and dielectric reciprocity
  before equivalent-atom transfer; full-set denoising also enforces equivalent-atom
  relationships. Residuals compare corrected back-contractions with measured seeds.
- Public broadening settings use `broadening_hwhm` in cm^-1, default 1.0.
- Removed symbolic tensor forms and polarization-angle scans from the public API
  and CLI. Fixed-polarization spectra and numerical symmetry diagnostics remain.
- Collective `mode-displacements` now defaults to explicit +/-Q pairs for
  every Raman-active Gamma mode, including every degenerate partner. Activity
  selection uses validated symmetry even without displacement reduction;
  inactive modes are no longer included by the unreduced automatic path.
  Partner/sign reduction requires `--reduce`. Unreduced mode-tensor extraction
  preserves the directly measured finite differences without orbit averaging.
- VASP dielectric processing now uses `vasprun.xml` exclusively. All nine
  named Cartesian components must be present; transpose elements are not
  synthesized. The former OUTCAR dielectric reader has been removed. Update
  old inputs to `dielectric_source: vasprun.xml`; six-component XML exports
  are intentionally rejected.
- Atomic derivatives always retain a `raw_directional_derivative_*.json` output. Site projection adds
  `corrected_` derivatives and diagnostics. Atomic spectrum generation now
  also preserves separately named raw and corrected tensors/intensities.
- Raman tensors have energy-specific filenames. The legacy `Raman_tensor`
  alias remains the last processed energy; unprefixed intensity files remain
  the selected raw or corrected result.
- Added experimental collective-mode generation, partner/sign reduction,
  central-difference reconstruction, and optional joint mode-covariance
  projection. Raw mode tensors remain available in NPZ format.
- Broadened spectra use linewidth-aware sampling (at most FWHM/20 spacing)
  instead of a fixed 2,000-point grid, and reject nonpositive linewidths.

These changes do not alter completed DFT outputs, launch simulations, or
demonstrate numerical convergence. Symmetry projection removes only the
symmetry-inconsistent component of error.
