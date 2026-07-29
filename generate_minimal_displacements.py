import numpy as np
import os
import shutil
import sys
from string import ascii_lowercase

from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp, write_vasp

from process_symmetry import read_contcar


def build_phonopy(contcar_path="CONTCAR", amplitude=0.03):
    """
    Builds a Phonopy object on the unit cell itself (supercell_matrix =
    identity), matching this pipeline's `phonopy: dim: "1 1 1"` convention --
    the Raman displacement problem is a single-cell, per-atom question, not a
    force-constant/supercell one, so no supercell expansion is needed.
    """
    cell = read_vasp(contcar_path)
    ph = Phonopy(cell, supercell_matrix=np.eye(3, dtype=int), primitive_matrix=np.eye(3))
    ph.generate_displacements(distance=amplitude, is_diagonal=False)
    return ph


def run_generate(contcar_path="CONTCAR", amplitude=0.03, template_dir=None):
    """
    Phonopy-driven replacement for raman_dis/raman_dis_nosym: generates only
    the symmetry-independent atoms x symmetry-independent Cartesian
    directions actually needed (via Phonopy.generate_displacements, i.e.
    phonopy's own get_least_displacements -- not a hand-rolled reimplementation).
    This is generally *fewer* displacements than raman_dis currently produces,
    since raman_dis only reduces atom count (via the space group) and not
    per-atom direction count (via site symmetry) -- see symmetry_operation_matrices
    for the atom-level reduction it already does, which this script reuses
    conceptually via phonopy's identical atom_mapping.

    Writes, per displacement:
      raman_poscar_<label><suffix>/POSCAR   -- ready for a VASP run
    plus a `displacements.dat` bookkeeping file in the same format
    calculate_dielectric_derivatives.py already expects (extended to support
    a variable number of displacements per atom, not a fixed 6).
    """
    if not os.path.exists(contcar_path):
        print(f"***** {contcar_path} not found *****")
        sys.exit(1)

    ph = build_phonopy(contcar_path, amplitude)
    positions_frac, total_atoms = read_contcar(contcar_path)

    # Atom symbols, in POSCAR order (needed for VASP-ready directories and
    # for the same "Elem<n>" per-atom labels create_displacements.py uses).
    with open(contcar_path) as f:
        lines = f.readlines()
    symbols_line = lines[5].split()
    counts_line = lines[6].split()
    try:
        counts = [int(x) for x in counts_line]
        symbols = symbols_line
    except ValueError:
        counts = [int(x) for x in symbols_line]
        symbols = [f"Elem{i+1}" for i in range(len(counts))]
    full_symbols = []
    for sym, c in zip(symbols, counts):
        full_symbols.extend([sym] * c)

    atom_label_counts = {s: 0 for s in symbols}
    per_atom_label = []
    for s in full_symbols:
        atom_label_counts[s] += 1
        per_atom_label.append(f"{s}{atom_label_counts[s]}")

    write_vasp("ref_poscar.vasp", ph.unitcell)

    suffix_counter = {}
    disp_records = []  # (label, atom_idx_1based, frac_disp_vector)
    lattice = ph.unitcell.cell

    for entry in ph.dataset["first_atoms"]:
        atom0 = entry["number"]           # 0-based
        cart_disp = entry["displacement"]  # Cartesian, Angstrom
        label = per_atom_label[atom0]

        suffix_counter.setdefault(label, 0)
        suffix = ascii_lowercase[suffix_counter[label]]
        suffix_counter[label] += 1

        dirname = f"ra_pos_atom{atom0 + 1}{suffix}"
        os.makedirs(dirname, exist_ok=True)

        displaced_positions = ph.unitcell.scaled_positions.copy()
        frac_disp = cart_disp @ np.linalg.inv(lattice)
        displaced_positions[atom0] += frac_disp

        displaced_cell = ph.unitcell.copy()
        displaced_cell.scaled_positions = displaced_positions
        write_vasp(os.path.join(dirname, "POSCAR"), displaced_cell)

        if template_dir:
            for fname in ("INCAR", "KPOINTS", "POTCAR"):
                src = os.path.join(template_dir, fname)
                if os.path.exists(src):
                    shutil.copyfile(src, os.path.join(dirname, fname))

        disp_records.append((label, atom0 + 1, frac_disp))
        print(f"  wrote {dirname}/POSCAR  (atom {label}, "
              f"cart disp {cart_disp.round(5).tolist()})")

    with open("displacements.dat", "w") as f:
        f.write("displacements for VASP. System (Phonopy site-symmetry-minimal set)\n")
        f.write(f"{len(disp_records):6d}   Number of displacements\n")
        for label, atom_idx, frac in disp_records:
            f.write(f"{label:<5s}{atom_idx:6d}{frac[0]:12.6f}{frac[1]:12.6f}{frac[2]:12.6f}\n")
        unique_atoms = sorted(set((r[0], r[1]) for r in disp_records), key=lambda x: x[1])
        f.write(f"{len(unique_atoms):6d}{total_atoms:6d}     Number of atoms in SC\n")
        for label, atom_idx in unique_atoms:
            f.write(f"{label:<5s}{atom_idx:6d}\n")

    n_full = total_atoms * 6
    print(f"\nGenerated {len(disp_records)} displacements for "
          f"{len(unique_atoms)} symmetry-independent atom(s) "
          f"(vs. {n_full} for a full, unreduced 6-per-atom set: "
          f"{100*len(disp_records)/n_full:.0f}% of the brute-force cost).")
    print("Wrote ref_poscar.vasp and displacements.dat.")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contcar", default="CONTCAR")
    ap.add_argument("--amplitude", type=float, default=0.03)
    ap.add_argument("--template-dir", default=None,
                     help="directory to copy INCAR/KPOINTS/POTCAR from into each raman_poscar_* dir")
    args = ap.parse_args()
    run_generate(args.contcar, args.amplitude, args.template_dir)
