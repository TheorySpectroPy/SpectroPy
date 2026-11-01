# Collective mode displacements and degeneracy symmetry

Experimental. No existing production inputs are changed.
No Fortran sources are modified. DFT is still run by the external workflow.

## Convention and Liangbo compatibility

Liangbo's raman_dis_mode.f90 reads CONTCAR and band.yaml, prompts for a
1-based mode index and Q (Angstrom sqrt(amu)), and writes:

- POSCAR_mode6a / POSCAR_mode6b: equilibrium Cartesian coordinates plus/minus
  Q times the mass-normalized phonon eigenvector e/sqrt(m).
- mode_index: three columns (mode index, frequency in cm^-1, Q).
- runRA: stages ra_POSCAR_mode* directories.

SpectroPy uses the SAME normal-coordinate convention, not a fixed maximum
atomic displacement. Atom registration is performed before displacement;
no additional mass factor belongs in the finite difference with denominator
2Q. JSON retains full Q and eigenvector precision, while mode_index keeps
Liangbo's human-readable field widths.

For an original-mode comparison:

    spectropy mode-displacements --modes 6 --q 0.01 --no-reduction

The generated directory contains the same POSCAR_mode6a/b interchange
filenames, mode_index, runRA, and prepared ra_POSCAR_mode6a/b/POSCAR.
SpectroPy writes the interchange POSCARs in Cartesian form like Liangbo.
Mode-feature validation uses VASP inputs; synthetic response tests are not
a substitute for independent DFT comparisons.

Verification against the unchanged, existing Fortran executable:
MoS2 mode 6 (E') and mode 8 (A1'), both Q=0.01, agreed within 4.6e-10
Angstrom in every Cartesian coordinate and exactly in the mode_index row.
The bound is consistent with the Fortran output's nine decimal places.
The comparison uses consistently ordered structure/phonon inputs.

## Reduced planning

    spectropy mode-displacements --q 0.01
    spectropy mode-displacements --modes 6 --q 0.01

The reduced path requires band.yaml, irreps.yaml and symmetry, plus
CONTCAR/ref_poscar.vasp for the VASP development path.
It calculates representation matrices in the actual registered eigenvector
basis, checks subspace closure and mass preservation, and groups modes from
Phonopy's irrep data, not from frequency proximity alone.
Each group is also required to be degenerate within 0.01 cm^-1.

Default automatic selection includes positive-frequency Raman-allowed
symmetric-tensor groups, excluding frequencies below 1 cm^-1. Activity is
determined from the rank of the tensor-intertwining equations, not irrep
name strings. Explicit --modes overrides activity/frequency filtering and
expands a requested mode to its entire partner group. In --no-reduction,
explicit original modes are kept separate.

For coefficients c in a partner basis, the collective pattern is
Q sum_p c_p e_p/sqrt(m). Its orbit coefficients are Gamma(g)c.
The planner chooses enough seed orbits to span the partner space. It tries
real -1 eigenvectors of representation matrices to reuse an opposite sign.
This is a deterministic greedy rank/cost search, NOT a proof of globally
minimal signed cost or optimal conditioning.

IMPORTANT: a reduced file named POSCAR_mode6a can represent a ROTATED
combination of modes 6 and 7. Its authoritative specification is the
mode_manifest.json seed coefficients, not the mode number alone.
The mode_index row labels the representative group; never feed it to an
older single-mode postprocessor while ignoring this metadata.

If g maps the complete collective pattern to its negative, the negative
dielectric response is S(g) epsilon(+Q) S(g)^T. Otherwise both signs are
computed explicitly. --keep-both-signs disables this reuse while retaining
partner reduction.

## Direct tensor reconstruction

    spectropy mode-tensors --source vasprun.xml

The manifest records the calculator and displacement basis.
For each seed, form [epsilon(+Q)-epsilon(-Q)]/(2Q). Rotate that response
through the same operations used for its partner coefficients, then solve
the full-rank partner system. Multiply by the cell volume/(4 pi) once.
This returns all requested partner tensors in the ORIGINAL Phonopy basis,
including their complex components. No atomic derivative reconstruction is
used. Only central differences are implemented for this feature.

The orbit fit enforces symmetry relations of the group. Its transformed
observations are not independent DFT measurements and do not imply a
statistical gain from independent noise averaging. Exact symmetry-sign reuse
does not make truncation errors inside the allowed subspace disappear.

Outputs in a separate mode_results directory:

- mode_tensors_<energy>eV.npz: exact mode IDs and complex tensor arrays.
- Raman_tensor_<energy>eV: mode index, THz, cm^-1 and complex matrix rows,
  following the existing Raman table layout but with 16-digit values.
- Raman_intensity_complex_<energy>eV: cm^-1, polarized contraction squared,
  label column; only requested/reconstructed modes, not a padded full spectrum.
- mode_covariance_residuals_<energy>eV.csv: absolute RMS mismatch per mode.

Dielectric tensors must contain every Cartesian component. Missing components
are never filled with zero. VASP uses the XML reader, preserving all nine
named components. Generated structures are checked
against the stored collective displacement before reading dielectric data.
Existing output directories are never silently overwritten.

## Optional mode-level denoising (VASP)

    spectropy mode-tensors --source vasprun.xml --denoise --output projected_mode_results

The complete partner set is projected jointly by the group average
Rbar_k = mean_g sum_p Gamma_kp(g) S_g R_p S_g^T.
The complete symmetry group and matching Gamma representation are validated
against the original band.yaml and current reference/symmetry/irreps.yaml.
Partial partner groups are rejected. Complex tensor values are preserved.
Do not average each partner independently as if it were totally symmetric.

Raw tensors are always saved in raw_mode_tensors_<E>eV.npz. With --denoise,
corrected_mode_tensors_<E>eV.npz and mode_projection_residuals_<E>eV.csv
are also written. The normal tensor/intensity outputs then use corrected
tensors. The CSV reports raw norm, removed norm, and absolute Neumann RMS
before/after correction, without intensity weighting. Vanishing corrected
residual is enforced consistency, not proof of physical accuracy. Reduced
orbit reconstruction can already enforce the same constraints, leaving
little further correction. Current intensity output requires specific
polarization; forward differences are not supported here.

## Workflow integration

Existing production simulations are unchanged.
VASP mode directories can be staged with the generated runRA and run by
the usual external VASP launcher; automatic VASP workflow array integration
is not yet implemented. Do not change a production material to this method
without choosing a fresh directory and verifying the required backend route.

## Scope

Symmetry operations must preserve the equilibrium ELECTRONIC state, including
magnetic order, and atomic masses. Symmetric dielectric tensors are assumed
for automatic Raman-activity selection. Complex irreps needing realification,
non-Gamma displacements, and general magnetic/corepresentation selection
are not supported. Labels or frequency degeneracy alone do not establish
valid partner relations.

Mode reduction is not a universal improvement over atomic SpectroPy
reduction; compare signed costs for the same structure and requested modes.
