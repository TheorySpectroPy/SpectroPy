# Change notes

## Unreleased — VASP workflow updates

- VASP dielectric processing now uses `vasprun.xml` exclusively. All nine
  named Cartesian components must be present; transpose elements are not
  synthesized. The former OUTCAR dielectric reader has been removed. Update
  old inputs to `dielectric_source: vasprun.xml`; six-component XML exports
  are intentionally rejected.
- Atomic derivatives always retain a `raw_` output. Site projection adds
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
