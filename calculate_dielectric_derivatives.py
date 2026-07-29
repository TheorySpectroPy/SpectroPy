import numpy as np
import os
import sys
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict

from process_symmetry import (
    read_contcar,
    read_symmetry_file,
    read_site_symmetries,
    compute_mapping_matrices,
)
from reconstruct_dielectric_derivatives import (
    frac_to_cart_rotation,
    compute_atom_tensor,
    expand_to_equivalent_atoms,
)


def read_diel_from_xml(xml_path, target_frequency):
    """
    Parses a vasprun.xml file to find the complex dielectric tensor
    at the frequency closest to the target frequency.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except (ET.ParseError, FileNotFoundError):
        print(f"Error: Could not parse or find {xml_path}")
        return None, None

    try:
        diel_func = root.find("./calculation/dielectricfunction")
        imag_part = diel_func.find('./imag/array/set')
        real_part = diel_func.find('./real/array/set')

        frequencies = []
        eps_imag_components = []
        eps_real_components = []

        for r in imag_part.findall('r'):
            data = list(map(float, r.text.split()))
            frequencies.append(data[0])
            eps_imag_components.append(data[1:])

        for r in real_part.findall('r'):
            data = list(map(float, r.text.split()))
            eps_real_components.append(data[1:])

        frequencies = np.array(frequencies)
        eps_imag_components = np.array(eps_imag_components)
        eps_real_components = np.array(eps_real_components)

        idx = np.argmin(np.abs(frequencies - target_frequency))

        ri = eps_real_components[idx]
        im = eps_imag_components[idx]

        eps_real = np.array([[ri[0], ri[3], ri[5]], [ri[3], ri[1], ri[4]], [ri[5], ri[4], ri[2]]])
        eps_imag = np.array([[im[0], im[3], im[5]], [im[3], im[1], im[4]], [im[5], im[4], im[2]]])

        return eps_real, eps_imag

    except AttributeError:
        print(f"Error: Could not find dielectric data in {xml_path}")
        return None, None


def run_generate_derivatives():
    """
    Collects VASP outputs and computes each atom's full dielectric-derivative
    tensor D_ijk, supporting a VARIABLE number of displacements per atom (not
    the fixed 6-per-atom the previous version of this script assumed) --
    generate_minimal_displacements.py's phonopy-driven set can have as few
    as 2 displacements for a highly-symmetric atom. Missing +/- pairs and
    entirely-skipped Cartesian axes are reconstructed via this atom's own
    site symmetry (compute_atom_tensor), and symmetry-equivalent atoms not
    displaced at all (e.g. MoS2's second S atom) are filled in via
    process_symmetry.py's atom-mapping matrices (expand_to_equivalent_atoms)
    -- this is the "full tensor reconstruction... not yet implemented" gap
    calculate_spectrum.py previously flagged.
    """
    required_files = ["ref_poscar.vasp", "displacements.dat", "CONTCAR", "symmetry"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    print("Reading displacements.dat to determine which files to collect...")
    with open("displacements.dat", 'r') as f:
        displacements_lines = f.readlines()
    num_disps = int(displacements_lines[1].split()[0])

    displacements = []
    for line in displacements_lines[2: 2 + num_disps]:
        parts = line.split()
        displacements.append({
            'label': parts[0],
            'index': int(parts[1]),
            'vector': np.array(list(map(float, parts[2:5]))),
        })

    print("Collecting vasprun.xml files...")
    os.makedirs("vasprun", exist_ok=True)

    alphabet = "abcdefghijklmnopqrstuvwxyz"
    disp_counter = defaultdict(int)
    suffixes = []
    for d in displacements:
        suffixes.append(alphabet[disp_counter[d['label']]])
        disp_counter[d['label']] += 1

    for disp_data, suffix in zip(displacements, suffixes):
        label = disp_data['label']
        source_dir_name = f"raman_poscar_{label}{suffix}"
        source_path = os.path.join(source_dir_name, "vasprun.xml")
        dest_filename = f"{label}{suffix}.xml"
        dest_path = os.path.join("vasprun", dest_filename)
        try:
            shutil.copyfile(source_path, dest_path)
        except FileNotFoundError:
            print(f"Error: Source file not found: {source_path}")
            print("Please ensure all VASP calculations have finished.")
            sys.exit(1)

    print("File collection complete. vasprun directory is ready.")
    print("-" * 40)

    print("Program <calculate_dielectric_derivatives.py>")
    print("Calculates derivatives of dielectric tensors from the collected vasprun.xml files.")
    try:
        target_frequency = float(input("Please input the laser frequency (eV): "))
    except ValueError:
        print("Invalid frequency. Exiting.")
        sys.exit(1)
    print("-" * 40)

    print("Reading ref_poscar.vasp...")
    with open("ref_poscar.vasp", 'r') as f:
        lines = f.readlines()
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])

    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[8:8 + total_atoms]])
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        total_atoms = sum(atom_counts)
        positions = np.array([list(map(float, l.split()[:3])) for l in lines[7:7 + total_atoms]])

    print("Processing vasprun.xml files from vasprun directory...")
    epsilon_log_filename = f"EPSILON_{target_frequency:.2f}.dat"
    per_atom_raw = defaultdict(list)  # atom_idx (1-based) -> [(cart_dir, eps_real, eps_imag), ...]

    with open(epsilon_log_filename, 'w') as f_log:
        f_log.write(f"Dielectric tensors for target frequency {target_frequency:.4f} eV\n")
        f_log.write("File                 Positions (frac)              Displacement (Angstrom)\n")
        f_log.write("-" * 80 + "\n")

        for disp_data, suffix in zip(displacements, suffixes):
            label = disp_data['label']
            xml_filename = f"{label}{suffix}.xml"
            xml_path = os.path.join("vasprun", xml_filename)

            eps_real, eps_imag = read_diel_from_xml(xml_path, target_frequency)
            if eps_real is None:
                sys.exit(1)

            cart_disp_vec = disp_data['vector'] @ lattice_vectors
            per_atom_raw[disp_data['index']].append((cart_disp_vec, eps_real, eps_imag))

            pos_str = " ".join([f"{p:8.5f}" for p in positions[disp_data['index'] - 1]])
            disp_str = " ".join([f"{d:8.5f}" for d in cart_disp_vec])
            f_log.write(f"{xml_filename:<12s}   {pos_str}      {disp_str}\n")
            for j in range(3):
                real_str = " ".join([f"{x:10.6f}" for x in eps_real[j]])
                imag_str = " ".join([f"{x:10.6f}" for x in eps_imag[j]])
                f_log.write(f"  real: {real_str} | imag: {imag_str}\n")
            f_log.write("\n")

    print(f"Finished processing XMLs. Log written to {epsilon_log_filename}")

    print("Reading symmetry data (site symmetries + atom mapping)...")
    site_syms = read_site_symmetries("symmetry")
    space_rotations, space_translations, atom_mapping = read_symmetry_file("symmetry")
    _, _, mapping_matrices_frac = compute_mapping_matrices(
        positions, space_rotations, space_translations, atom_mapping
    )
    mapping_matrices_cart = {
        pair: frac_to_cart_rotation(R, lattice_vectors)
        for pair, R in mapping_matrices_frac.items()
    }

    print("Reconstructing full per-atom D_ijk tensors from the displacement set actually run "
          "(filling in any skipped +/- pairs or axes via this atom's own site symmetry)...")
    D_real_by_independent = {}
    D_imag_by_independent = {}
    for atom_idx, raw in per_atom_raw.items():
        ss_frac = site_syms[atom_idx]["rotations"]
        ss_cart = np.array([frac_to_cart_rotation(R, lattice_vectors) for R in ss_frac])
        measured_real = [(d, er) for d, er, ei in raw]
        measured_imag = [(d, ei) for d, er, ei in raw]
        D_real_by_independent[atom_idx] = compute_atom_tensor(measured_real, ss_cart)
        D_imag_by_independent[atom_idx] = compute_atom_tensor(measured_imag, ss_cart)

    print("Expanding to symmetry-equivalent atoms not directly displaced...")
    D_real_all = expand_to_equivalent_atoms(D_real_by_independent, atom_mapping, mapping_matrices_cart)
    D_imag_all = expand_to_equivalent_atoms(D_imag_by_independent, atom_mapping, mapping_matrices_cart)

    dielectric_derivatives_filename = f"dielectric_derivatives_{target_frequency:.2f}"
    print(f"Writing final derivatives of dielectric tensors to {dielectric_derivatives_filename}...")

    pi = np.pi
    cell_volume = np.linalg.det(lattice_vectors)
    conv_factor = cell_volume / (4 * pi)

    all_atom_indices = sorted(D_real_all.keys())
    with open(dielectric_derivatives_filename, 'w') as f:
        f.write(f"! derivatives of dielectric tensors calculated for laser frequency {target_frequency:.4f} eV\n")
        f.write("! from VASP vasprun.xml files\n")
        f.write("!\n! Unit cell matrix:\n")
        for vec in lattice_vectors.T:
            f.write(f"!   {vec[0]:21.16f} {vec[1]:21.16f} {vec[2]:21.16f}\n")
        f.write("!\n! derivatives of dielectric tensors\n")

        f.write("! The Real Part of dielectric tensor derivatives:\n")
        for atom_idx in all_atom_indices:
            pos = positions[atom_idx - 1]
            f.write(f"      atom{atom_idx} {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}\n")
            D = D_real_all[atom_idx] * conv_factor
            tensor_to_print = np.hstack((D[:, :, 0], D[:, :, 1], D[:, :, 2]))
            for row in tensor_to_print:
                f.write("".join([f"{x:16.4f}" for x in row]) + "\n")

        f.write("\n! The Imaginary Part of dielectric tensor derivatives:\n")
        for atom_idx in all_atom_indices:
            pos = positions[atom_idx - 1]
            f.write(f"      atom{atom_idx} {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}\n")
            D = D_imag_all[atom_idx] * conv_factor
            tensor_to_print = np.hstack((D[:, :, 0], D[:, :, 1], D[:, :, 2]))
            for row in tensor_to_print:
                f.write("".join([f"{x:16.4f}" for x in row]) + "\n")

    print("\ncalculate_dielectric_derivatives.py finished successfully.")
    print(f"Wrote data for {len(all_atom_indices)} atom(s) "
          f"({len(per_atom_raw)} directly computed, "
          f"{len(all_atom_indices) - len(per_atom_raw)} filled in by symmetry).")


if __name__ == "__main__":
    run_generate_derivatives()
