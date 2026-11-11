# MoS2 collective mode example

Provide consistently ordered CONTCAR (VASP), plus Phonopy
band.yaml, irreps.yaml and symmetry. No electronic results are bundled.
Use the neighboring input file for energies and xx polarization.

    spectropy mode-displacements --q 0.01

The default enumerates every positive-frequency Raman-active mode above
1 cm^-1, every degenerate partner separately, with both signs. It does not
reuse signs or reconstruct partners. Q=0.01 means Angstrom sqrt(amu), not a
fixed 0.01 Angstrom atomic displacement. For the local validation fixture this selects
modes 4, 5, 6, 7, 8 (E'', E', A1'), including out-of-plane optical channels
even though the example input requests xx polarization.

To select an E' pair explicitly (verify its mode numbers in your own dataset):

    spectropy mode-displacements --modes 6 7 --q 0.01 --output explicit_pair

Both modes are directly displaced. Explicit selection overrides Raman-activity
and frequency filtering. Every selected mode requires both signs.

The POSCAR_mode6a/b files and mode_index columns follow Liangbo's
raman_dis_mode conventions. No Fortran code is called by SpectroPy.

Stage your optical INCAR, KPOINTS and POTCAR from the calculation directory:

    bash mode_displacements/runRA

Run VASP externally in each mode_displacements/ra_POSCAR_mode*/ directory.
Do not relax displaced structures. SpectroPy does not submit jobs or choose
electronic convergence settings.

After external optical calculations finish in the generated directories:

    spectropy mode-tensors --directory mode_displacements --output mode_results

vasprun.xml is the only VASP dielectric source; all nine components are retained.
Each output tensor comes directly from its calculated positive/negative pair.
See [the main README](../../README.md#mode-displacements) for inputs and units.
