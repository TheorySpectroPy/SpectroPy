# MoS₂: from atomic displacements to a Raman spectrum

This tutorial starts from an equilibrium structure and precomputed Γ-point
phonons, **before calculating dielectric derivatives**. You will generate
displaced structures, supply their dielectric responses, reconstruct atomic
derivatives, and calculate a polarized Raman spectrum.

**Data pending:** the scientific input/output files have not yet been selected
or bundled. The `data/` directories contain explicit placeholders, not usable
calculation data. Installation and help commands work now; stop before step 2
until the starting files are supplied. The dielectric-output bundle required
for the no-DFT route will be supplied separately. No numerical result or
convergence claim is attached to the current template.

## Scope: independent of your execution workflow

SpectroPy does not choose a functional, potentials, k-points, executable,
cluster scheduler, or relaxation/phonon workflow. This example uses its
**VASP file interface**, so it is not file-format agnostic: the structural
input is `CONTCAR`, and dielectric data come from nine-component `vasprun.xml`
files. No Fortran program or production workflow generator is needed here.
DFT execution is an explicit external step, with a prerecorded-data alternative.

Relaxation and harmonic-force calculations are prerequisites, not performed
by this tutorial. The structure, phonons, symmetry data, and optical calculations
must describe the same equilibrium system and compatible Cartesian frame.

## 1. Install and choose a clean working directory

In a terminal at the **SpectroPy repository root**, run:

```bash
python -m pip install -e '.[phonopy,plot]'
spectropy --help
spectropy displacements --help

# Keep source/reference files separate from generated results.
TUTORIAL_DIR="$PWD/examples/MoS2_From_Displacements"
RUN_DIR=$(mktemp -d "${TMPDIR:-/tmp}/spectropy-mos2-XXXXXX")
cp "$TUTORIAL_DIR/input" "$RUN_DIR/input"
cd "$RUN_DIR"
```

Keep this terminal open so the two path variables remain defined. `RUN_DIR`
is a new temporary directory; use a new durable directory instead if your
machine periodically purges temporary files. No command below deletes data.

## 2. Supply and inspect the starting files — placeholder

We will provide these files in `data/starting/`:

| File | Purpose |
| --- | --- |
| `CONTCAR` | Relaxed three-atom monolayer MoS₂ reference cell. |
| `band.yaml` | Phonopy Γ-point frequencies, masses, structure, and eigenvectors. |
| `symmetry` | Phonopy YAML containing space-group operations, atom mapping, and site symmetries. |
| `irreps.yaml` | Mode labels and partner groups for the optional mode-covariance check. |

The exact calculation provenance and versions will accompany the data. Do not
substitute the older tutorial's derivatives: generating derivatives is the
point of this example. Do not combine unrelated structure and phonon datasets.

Once the files exist, copy them with this guarded command:

```bash
(
  set -eu
  for name in CONTCAR band.yaml symmetry irreps.yaml; do
    test -s "$TUTORIAL_DIR/data/starting/$name" || {
      echo "Missing starting file: $name — supply the tutorial data first." >&2
      exit 1
    }
  done
  for name in CONTCAR band.yaml symmetry irreps.yaml; do
    cp "$TUTORIAL_DIR/data/starting/$name" .
  done
)
```

**Stop if this command reports a missing file.** The placeholders cannot be
renamed into scientific inputs. Inspect `input` before proceeding: it selects
xx polarization, a 0.03 Å displacement length, central differences, site
projection, and a static-limit dielectric response at 0.00 eV. These are
tutorial choices, not a prescription for converged calculations. No `cutoff`
is set, so intensities are reduced squared Raman contractions without the
optional Bose/phonon-frequency factor.

## 3. Generate symmetry-reduced atomic displacements

Run from `RUN_DIR`:

```bash
spectropy displacements --mode fractional
spectropy symmetry
```

The fractional direction basis is converted into physical displacements of
the requested length. The output includes:

```text
ref_poscar.vasp
displacements.dat
atomic_displacement
fractional_basis.json
symmetry_operation_matrices
displacements/
  pos_atom...                 # generated structure files
  ra_pos_atom.../
    POSCAR                   # structure to use for that optical calculation
```

The ellipses illustrate naming patterns, not literal names. Read the actual
names and displacements in `displacements.dat`; do not hard-code their count.
Equivalent atoms/directions and, where justified, opposite signs are supplied
by symmetry. A generated directory does not necessarily correspond to an
independent displacement direction. Do not change these POSCARs or regenerate
with different settings after calculating their dielectric responses.

## 4. Supply dielectric responses: choose one route

### A. Run your own optical calculations

For each `displacements/ra_pos_atom*/` directory:

1. Keep its generated `POSCAR` unchanged. This is a frozen displacement, not a
   new structure to relax.
2. Supply your independently chosen VASP optical inputs and licensed potentials.
   SpectroPy does not create `INCAR`, `KPOINTS`, `POTCAR`, or submission scripts.
   Use consistent electronic settings and potentials across the displaced set.
3. Run VASP through your usual local or cluster launcher. No specific `srun`,
   `mpirun`, executable path, or resource allocation is prescribed here.
4. Verify electronic convergence and successful completion. Confirm that
   `vasprun.xml` contains the density-density dielectric function, with real
   and imaginary data covering each requested energy.

**All nine named components are required:** `xx`, `xy`, `xz`, `yx`, `yy`,
`yz`, `zx`, `zy`, and `zz`, plus `energy`. Not every VASP build exports this
layout. SpectroPy rejects six-component tables; it does not invent transpose
entries or fall back to OUTCAR. Check one completed calculation before running
the entire set. A parseable XML file alone does not establish SCF convergence.

### B. Use the prerecorded outputs — placeholder

The future `data/dielectrics/` bundle will contain:

```text
displacements.dat
ra_pos_atom.../
  POSCAR
  vasprun.xml
```

It must match the generated displacement set exactly. After the bundle is
supplied, the following performs a strict metadata/structure comparison and
then copies only the XML files. It intentionally stops on any mismatch or
pre-existing destination XML; do not bypass a mismatch by renaming files.

```bash
(
  set -eu
  bundle="$TUTORIAL_DIR/data/dielectrics"
  test -s "$bundle/displacements.dat"
  cmp displacements.dat "$bundle/displacements.dat"
  for directory in displacements/ra_pos_atom*; do
    test -d "$directory"
    name=$(basename "$directory")
    cmp "$directory/POSCAR" "$bundle/$name/POSCAR"
    test -s "$bundle/$name/vasprun.xml"
    test ! -e "$directory/vasprun.xml"
  done
  for directory in displacements/ra_pos_atom*; do
    name=$(basename "$directory")
    cp "$bundle/$name/vasprun.xml" "$directory/vasprun.xml"
  done
)
```

A format-only mismatch can also fail `cmp`; ask for a compatible bundle rather
than assuming equivalence. The release bundle must be generated with the
documented SpectroPy version and the provided `input`. It is not available yet.

## 5. Reconstruct derivatives, keeping raw and denoised data

```bash
spectropy derivatives
```

For the supplied `input`, expect:

- `dielectric_tensor_0.00`: sampled dielectric responses.
- `raw_epsilon_derivative_0.00`: reconstructed derivatives before site averaging.
- `corrected_epsilon_derivative_0.00`: site-projected derivatives.
- `epsilon_derivative_0.00`: the active corrected field for downstream consumers.
- `site_projection_residuals_0.00.csv`: per-atom raw-to-corrected differences.

“Raw” still includes symmetry reconstruction. Site averaging acts on the full
three-slice per-atom tensor, not just forbidden entries. Projection enforces
symmetry consistency; it cannot detect or remove every numerical error. The
derivative files include the cell-volume/(4π) factor; do not multiply it again.

The default path reconstructs direction bases explicitly. The optional
`spectropy derivatives --saved-basis` validates and reuses `fractional_basis.json`.
Changing the displacement length requires new dielectric calculations, not
just editing metadata.

## 6. Calculate and plot the xx spectrum

```bash
spectropy spectrum
spectropy plot
```

Tensor and intensity tables are saved with `raw_` and `corrected_` prefixes.
Unprefixed tables select the corrected result because projection is enabled.
The selected outputs include:

- `Raman_tensor_0.00eV` (and the legacy `Raman_tensor` alias).
- `Raman_intensity_complex_0.00eV`.
- `Raman_intensity_complex_broadening_0.00eV`.
- `Raman_plot_styled_Raman_intensity_complex_0.00eV.png`.

The standard plot command plots the **selected**, unprefixed data, not both
raw and corrected variants. The PNG is normalized for display; the broadened
numerical spectrum is unnormalized. A 1 cm⁻¹ Lorentzian FWHM controls display
broadening, not the accuracy of the electronic calculation. No reference peak
positions or expected intensity ratios are supplied until the data are finalized.

## 7. Optional checks and extensions

```bash
spectropy check --kind mode --plot
```

This compares mode tensors against their symmetry transformation law using
the phonon partner groups. A small residual is a consistency check, not proof
of physical accuracy. Read the generated diagnostics before interpreting it.

To examine xy polarization, change the second line of `input` to `0.0 1.0 0.0`,
then rerun `spectropy spectrum` and `spectropy plot`. This needs no new DFT
calculations, but overwrites the same spectrum filenames: preserve the xx
outputs first, or use a separate working copy. Likewise, additional photon
energies need no new DFT calculations only if the existing dielectric tables
cover them. Rerun derivatives, spectrum, and plot after changing energies.

For more detail, see [site projection](../../docs/site_projection.md) and
the [main workflow documentation](../../README.md). Collective-mode displacements
are a different route; see the [separate example](../MoS2_Mode_Displacements/README.md).
