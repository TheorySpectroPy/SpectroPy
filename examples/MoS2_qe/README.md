# MoS₂: Quantum ESPRESSO — WIP

**Work in progress: documentation starter only, not a runnable example.**
The supported SpectroPy commands currently use the VASP file interface.
This directory does not enable a QE backend or introduce a working `--qe` flag.
For the complete, tested example, use [MoS2_vasp](../MoS2_vasp/README.md).

## Planned example

Follow the same sequence as the VASP tutorial: generate symmetry-reduced
atomic displacements, run external electronic calculations, reconstruct and
denoise the atomic derivatives, and calculate Raman tensors and spectra.

| Required input/output | Status |
| --- | --- |
| `qe.in` with the equilibrium structure and electronic settings | To be supplied and validated. |
| Consistently generated Γ-point phonons and site symmetries | To be supplied. |
| Displaced QE input files | Awaiting interface integration. |
| `epsilon.in` and matching dielectric-response outputs | To be supplied and validated. |
| Pseudopotential provenance and convergence settings | To be documented. |
| Reproducible Raman tensors and a spectrum plot | Pending an end-to-end test. |

Validation must check coordinate units and atom ordering, dielectric component
conventions, displacement lengths, and the full-versus-minimal comparison.
Do not substitute VASP benchmark XMLs or assume its settings transfer unchanged.
No QE numerical results, pseudopotentials, or experimental backend code are
distributed in this starter.
