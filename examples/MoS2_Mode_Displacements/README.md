# MoS2 collective mode example

Provide consistently ordered CONTCAR (VASP), plus Phonopy
band.yaml, irreps.yaml and symmetry. No electronic results are bundled.
Use the neighboring input file for energies and xx polarization.

    spectropy mode-displacements --q 0.01

To select an E' pair (verify its mode numbers in your own phonon dataset):

    spectropy mode-displacements --modes 6 --q 0.01

Mode 7 is then reconstructed jointly with mode 6, not dropped.
Mode filenames can refer to rotated partner combinations; inspect
mode_displacements/mode_manifest.json before running any DFT.

For the Fortran-compatible explicit original-mode stencil:

    spectropy mode-displacements --modes 6 --q 0.01 --no-reduction --output unreduced

The POSCAR_mode6a/b files and mode_index columns follow Liangbo's
raman_dis_mode conventions. No Fortran code is called by SpectroPy.

After external optical calculations finish in the generated directories:

    spectropy mode-tensors --directory mode_displacements --output mode_results

vasprun.xml is the only VASP dielectric source; all nine components are retained. Add --denoise to
project complete partner groups jointly and retain raw/corrected tensors.
QE workflow development remains separate from this VASP mode example.
See ../../docs/mode_displacements.md for assumptions and limitations.
