# MoS₂

A complete analysis example using the **minimal Cartesian, 0.01 Å PBEsol
benchmark**. The three calculated structures and their full `vasprun.xml`
outputs are included. You do not need VASP to reproduce the tensor and spectrum
analysis below.

## Run the example

From the repository root:

```bash
python -m pip install '.[all]'
cd examples/MoS2
spectropy tensors
spectropy spectrum
spectropy plot
```

| Command | Output |
| --- | --- |
| `tensors` | `Raman_tensor_<E>eV`: final complex Cartesian mode tensors. |
| `spectrum` | `Raman_intensity_complex_<E>eV` and `Raman_intensity_polarization_averaged_<E>eV`. |
| `plot` | PNG plots and unnormalized two-column broadened spectra. |

`<E>` is `0.00`, `1.96`, or `2.33`. Settings are grouped in `input.yaml`:
central differences, symmetry denoising, xx polarization plus the in-plane
polarization-channel sum, and Lorentzian HWHM 1 cm⁻¹. `spectrum.cutoff: 0`
retains the benchmark's Bose/phonon-frequency factor for positive-frequency
modes. PNGs are normalized for display; numerical spectra are not.

The tensor stage reconstructs atomic derivatives, denoises them, and contracts
them with the registered phonon eigenvectors in memory. Per-atom and raw
directional debugging files are optional; see [Advanced Usage](../../docs.md).

## Included files

```text
MoS2/
├── input.yaml
├── CONTCAR
├── ref_poscar.vasp
├── INCAR
├── KPOINTS
├── band.yaml
├── irreps.yaml
├── symmetry
├── displacements.dat
└── displacements/
    ├── ra_pos_atom1a/
    │   ├── POSCAR
    │   ├── vasprun.xml
    │   ├── INCAR   -> ../../INCAR
    │   └── KPOINTS -> ../../KPOINTS
    ├── ra_pos_atom2a/          # Same four files
    └── ra_pos_atom2b/          # Same four files
```

| File | Contents |
| --- | --- |
| `CONTCAR`, `ref_poscar.vasp` | Equilibrium geometry and displacement reference. |
| `band.yaml`, `irreps.yaml` | Γ-point phonon eigenvectors/frequencies and mode labels. |
| `symmetry` | Phonopy atom mappings and site-symmetry operations. |
| `displacements.dat` | Actual minimal displacement plan; vectors are stored fractionally. |
| `ra_pos_atom*/POSCAR` | Three directly calculated, mixed-direction displaced structures. |
| `ra_pos_atom*/vasprun.xml` | Unmodified electronic outputs, including all nine dielectric components. |
| `INCAR`, `KPOINTS` | The benchmark's optical-calculation inputs, unchanged. |

Three calculation folders do **not** mean three Cartesian derivative slices:
the Mo atom has one calculated sign, with its opposite supplied by site symmetry;
the representative S atom has both signs. The remaining directions and equivalent
S atom are reconstructed. The bundled metadata, POSCARs, and XML files belong to
the same run; do not combine them with the former tutorial dataset.

## Provenance and calculation settings

Source dataset: `MoS2_benchmark/MoS2_PBESol_minimal_Cartesian_d0.01/raman`,
completed 5 November 2026. Structures, phonons, symmetry, displacement metadata,
and XML files are copied without alteration. `input.yaml` expresses the original
analysis choices using the current grouped schema. No calculation is rerun by
this example, and these settings are not a general convergence recommendation.

| Optical setting | Value |
| --- | --- |
| VASP | 6.4.3, build 26 August 2026 |
| Functional | PBEsol (`GGA=PS`) |
| Potentials | PAW_PBE Mo_pv 04Feb2005 (12 valence electrons), S 06Sep2000 (6 each) |
| `ENCUT`, `NBANDS` | 400 eV, 64 |
| K-point mesh | Γ-centered 24×24×1 |
| `EDIFF`, `ISMEAR`, `SIGMA` | 10⁻⁸ eV, 0, 0.01 eV |
| `ISTART`, `ICHARG`, `ISYM` | 0, 2, −1 |
| `NSW`, `LOPTICS` | 0, true |
| `NEDOS`, `OMEGAMAX` | 50001 requested, 50 eV |

## Run new electronic calculations instead

The included XML files already provide the electronic results. To create a new
calculation, copy only the small inputs into a **new** directory, then run
`spectropy displacements`. Supply your own licensed `POTCAR` first; generation
symlinks supplied `INCAR`, `KPOINTS`, and `POTCAR` into each calculation directory.
Verify the generated plan and run VASP externally without relaxing displaced
structures. Do not overwrite the bundled benchmark XMLs.

`POTCAR`, VASP binaries, wavefunctions, and cluster-specific submission scripts
are intentionally not distributed. Execution and convergence checks remain the
user's responsibility. See the [main workflow](../../README.md#workflow--usage-guide)
for a from-scratch demonstration.
