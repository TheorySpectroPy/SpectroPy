# Placeholder: raw displaced-cell dielectric outputs

No electronic outputs are supplied yet. Provide the exact displacement
metadata `displacements.dat` and one directory per generated ra_pos_atom*:

    ra_pos_atom.../POSCAR
    ra_pos_atom.../vasprun.xml

Generate this bundle from the tutorial's finalized starting data and input.
Keep the original printed components in both real and imaginary density-density
tables and include the requested 0.00 eV point. Nine components allow independent
transpose comparisons; six-component XML is accepted with a reciprocity warning.
Do not replace XML outputs with precomputed derivatives or manually fill missing
components in the source files. Do not redistribute licensed POTCAR files.

Record the executable version/build, optical settings, convergence/completion
evidence, energy grid, output precision, units, source permission, and checksums.
If output files are hosted externally, add the download URL, archive layout,
size, and checksum here. Retain matching generated POSCARs for verification.
The README's copy command expects files unpacked directly into this directory.
