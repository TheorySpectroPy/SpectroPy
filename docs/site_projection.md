# Representative-first site-symmetry projection

Optional SpectroPy input setting: `site_symmetry_projection: true`.
Default: false, preserving the uncorrected workflow.

For a newly generated production workflow, use the YAML boolean
`spectropy.site_symmetry_projection: true` in the material settings (nested
under `spectropy:`), then run `raman-generate`. Existing generated workflows
are not changed automatically. Preserve their raw outputs before opting in.

For a symmetry-reduced set, reconstruct each measured representative's complete
Cartesian rank-three dielectric derivative with the existing direction solver.
Average its tensor transforms over its site subgroup, then generate equivalent
atoms by the existing space-group transfer. Real and imaginary parts are treated
separately. Both central and forward differences are supported.

The average enforces all site-component relations, not just forbidden zeros.
It neither enforces dielectric-index reciprocity nor an acoustic sum rule.
Operations must preserve the electronic state; antiunitary operations are not
supported. Extra independently measured atoms are rejected rather than discarded.
Full/no-symmetry displacement sets are not supported by this option.

Each derivative run writes `raw_epsilon_derivative_<energy>` with the
unprojected full atomic field, and `epsilon_derivative_<energy>` with the field
used for spectra. When projection is enabled, it also writes
`corrected_epsilon_derivative_<energy>` and
`site_projection_residuals_<energy>.csv`. The latter lists each atom's
species, representative status, Frobenius norm of its complex raw minus
corrected rank-three tensor, and that norm multiplied by the unit-cell
volume divided by 4π (the scale in the derivative table files). The console
prints each residual and the unweighted mean over measured representatives.

`spectropy spectrum` always writes `raw_Raman_tensor_<energy>eV` and
`raw_Raman_intensity_complex_<energy>eV`. With projection enabled it additionally
writes the corresponding `corrected_` files. Polarization-averaged intensity
files, when requested, follow the same prefixes. Unprefixed spectrum files
remain the active selection (corrected when enabled, raw otherwise), and
`Raman_tensor` remains the legacy alias for the last processed energy.
Energy-specific tensor files prevent one energy from overwriting another.
If projection is enabled, missing raw or corrected derivative files are an
error: rerun the derivative stage rather than relabeling corrected data as raw.
Raw here means before site averaging; symmetry reconstruction is still used.
Previously generated files are not deleted when a setting is disabled.

The correction can
lower the Neumann residual by construction; this is not independent validation
of physical accuracy. Errors already satisfying symmetry survive.
