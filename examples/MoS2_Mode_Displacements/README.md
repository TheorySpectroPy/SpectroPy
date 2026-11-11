# MoS2 collective mode example

Provide consistently ordered CONTCAR (VASP), plus Phonopy
band.yaml, irreps.yaml and symmetry. No electronic results are bundled.
Use the neighboring `input.yaml` for mode selection, amplitude, energies, output
directories, and xx polarization.
Place your optical INCAR, KPOINTS, and POTCAR in the working directory before
generation; SpectroPy automatically symlinks them into each displacement directory.

    spectropy mode-displacements

The default enumerates every positive-frequency Raman-active mode above
1 cm^-1, every degenerate partner separately, with both signs. It does not
reuse signs or reconstruct partners. Q=0.01 means Angstrom sqrt(amu), not a
fixed 0.01 Angstrom atomic displacement. For the local validation fixture this selects
modes 4, 5, 6, 7, 8 (E'', E', A1'), including out-of-plane optical channels
even though the example input requests xx polarization.

To select an E' pair explicitly (verify its mode numbers in your own dataset):

Set `mode_indices: [6, 7]` and `mode_displacements_directory: explicit_pair` in
`input.yaml`, then run `spectropy mode-displacements` in a clean working directory.

Both modes are directly displaced. Explicit selection overrides Raman-activity
and frequency filtering. Every selected mode requires both signs.

The POSCAR_mode6a/b files and mode_index columns follow Liangbo's
raman_dis_mode conventions. No Fortran code is called by SpectroPy.

Verify the generated input symlinks. Missing inputs are reported; no separate
runRA staging script is generated or needed.
Run VASP externally in each mode_displacements/ra_POSCAR_mode*/ directory.
Do not relax displaced structures. SpectroPy does not submit jobs or choose
electronic convergence settings.

After external optical calculations finish in the generated directories:

    spectropy mode-tensors

vasprun.xml is the only VASP dielectric source; all printed components are retained.
Six-component files are accepted with a warning that missing transpose entries
are filled by reciprocity, rather than independently calculated.
Each output tensor comes directly from its calculated positive/negative pair.
See [advanced documentation](../../docs.md#collective-mode-displacements) for inputs and units.
