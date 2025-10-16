import numpy as np
import os
import sys
import yaml

def read_contcar(filepath="CONTCAR"):
    """Reads a VASP CONTCAR file to get structure information."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    comment = lines[0].strip()
    # Assuming lattice constant is 1.0, as vectors are usually scaled in CONTCAR
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    
    try:
        atom_counts = np.array(lines[6].split(), dtype=int)
        atom_symbols = lines[5].split()
        header = lines[:8]
        pos_start_idx = 8
    except (ValueError, IndexError):
        atom_counts = np.array(lines[5].split(), dtype=int)
        atom_symbols = [f"El{i+1}" for i in range(len(atom_counts))]
        header = lines[:7]
        pos_start_idx = 7
    
    total_atoms = sum(atom_counts)
    positions_frac = np.array([list(map(float, l.split()[:3])) for l in lines[pos_start_idx:pos_start_idx+total_atoms]])
    positions_cart = positions_frac @ lattice_vectors

    return {
        "lattice_vectors": lattice_vectors,
        "atom_symbols": atom_symbols,
        "atom_counts": atom_counts,
        "total_atoms": total_atoms,
        "positions_cart": positions_cart,
        "header_lines": header # For recreating supercell POSCAR
    }

def read_band_yaml(filepath="band.yaml"):
    """Parses a Phonopy band.yaml file to get phonon modes."""
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    phonons = data['phonon'][0]['band']
    n_modes = len(phonons)
    n_atoms = data['natom']
    masses = np.array([atom['mass'] for atom in data['atoms']])
    
    frequencies = np.array([p['frequency'] for p in phonons])
    eigenvectors_raw = np.array([p['eigenvector'] for p in phonons])
    
    # Eigenvectors are complex, take the real part and reshape
    eigenvectors = eigenvectors_raw[..., 0].reshape(n_modes, n_modes)
    
    # Normalize by mass to get real-space eigendisplacements
    masses_expanded = np.repeat(masses, 3)
    eigendisplacements = eigenvectors / np.sqrt(masses_expanded)
    eigendisplacements = eigendisplacements.reshape(n_modes, n_atoms, 3)
    
    return frequencies, eigendisplacements, masses

def write_vmd_script(filename, positions, displacements, l_cylinder, l_cone):
    """Writes the VMD drawing commands for a single mode."""
    with open(filename, 'w') as f:
        f.write("draw color blue\n")
        for i in range(len(positions)):
            pos = positions[i]
            disp = displacements[i]
            
            # Start and end points of the cylinder (arrow body)
            p1 = pos
            p2 = pos + l_cylinder * disp
            
            # End points of the cone (arrow head)
            p3 = p2 + l_cone * disp
            
            f.write(f"draw cylinder {{{p1[0]:.6f} {p1[1]:.6f} {p1[2]:.6f}}} {{{p2[0]:.6f} {p2[1]:.6f} {p2[2]:.6f}}} radius 0.1 resolution 30 filled yes\n")
            f.write(f"draw cone {{{p2[0]:.6f} {p2[1]:.6f} {p2[2]:.6f}}} {{{p3[0]:.6f} {p3[1]:.6f} {p3[2]:.6f}}} radius 0.2 resolution 30\n")

def run_visualization():
    """Main function to generate visualization scripts for phonon modes."""
    # --- Check for required input files ---
    required_files = ["CONTCAR", "band.yaml"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    # --- Read input files ---
    print("Reading CONTCAR and band.yaml...")
    structure = read_contcar()
    frequencies, eigendisps, masses = read_band_yaml()
    n_modes = len(frequencies)

    # --- User Prompts ---
    try:
        factor = float(input("Please input the normalized factor to get proper arrow length: "))
    except ValueError:
        print("Invalid factor. Exiting.")
        sys.exit(1)
        
    # Scale factor 
    factor *= np.sqrt(np.max(masses))
    l_cylinder = 4.0 * factor
    l_cone = 1.5 * factor

    # --- Generate VMD scripts for primitive cell ---
    print(f"Generating {n_modes} VMD scripts for primitive cell (mode*.vmd)...")
    for i in range(n_modes):
        filename = f"mode{i+1}.vmd"
        write_vmd_script(filename, structure["positions_cart"], eigendisps[i], l_cylinder, l_cone)

    # --- Write frequency list file ---
    if not os.path.exists("all_modes.txt"):
        print("Writing all_modes.txt...")
        with open("all_modes.txt", "w") as f:
            f.write("# unit: cm-1\n")
            thz_to_cm1 = 33.35641
            for i in range(n_modes):
                f.write(f"{i+1:4d}   {frequencies[i] * thz_to_cm1:10.4f}\n")

    # --- Optional: Supercell Visualization ---
    supercell_choice = input("Visualize in a supercell? (yes/no): ").lower()
    if supercell_choice == 'yes':
        try:
            supercell_dim_str = input("Please input the dimension of supercell (e.g., 2 2 1): ")
            supercell_dim = np.array(supercell_dim_str.split(), dtype=int)
        except (ValueError, IndexError):
            print("Invalid supercell dimension. Exiting.")
            sys.exit(1)

        print("Generating supercell and VMD scripts...")
        
        # Create supercell positions
        supercell_positions = []
        prim_pos = structure["positions_cart"]
        lv = structure["lattice_vectors"]
        
        for i in range(supercell_dim[0]):
            for j in range(supercell_dim[1]):
                for k in range(supercell_dim[2]):
                    translation = i * lv[0] + j * lv[1] + k * lv[2]
                    supercell_positions.append(prim_pos + translation)
        
        supercell_positions = np.concatenate(supercell_positions, axis=0)
        
        # Write POSCAR for the supercell
        with open("POSCAR_supercell", "w") as f:
            f.write("Supercell for visualization\n")
            f.write("1.0\n")
            supercell_lvs = structure["lattice_vectors"] * supercell_dim[:, np.newaxis]
            for vec in supercell_lvs:
                f.write(f"   {vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}\n")
            f.write("  " + " ".join(structure["atom_symbols"]) + "\n")
            supercell_counts = structure["atom_counts"] * np.prod(supercell_dim)
            f.write("  " + " ".join(map(str, supercell_counts)) + "\n")
            f.write("Cartesian\n")
            for pos in supercell_positions:
                f.write(f"   {pos[0]:.9f} {pos[1]:.9f} {pos[2]:.9f}\n")
                
        # Generate VMD scripts for supercell
        for i in range(n_modes):
            filename = f"mode_super{i+1}.vmd"
            # Replicate the primitive cell displacements for each translated image
            supercell_disps = np.tile(eigendisps[i], (np.prod(supercell_dim), 1))
            write_vmd_script(filename, supercell_positions, supercell_disps, l_cylinder, l_cone)
            
    print("\nvisualize_modes.py finished successfully.")

if __name__ == "__main__":
    run_visualization()