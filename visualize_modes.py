import numpy as np
import os
import sys
import yaml
import re

def read_contcar(filepath="CONTCAR"):
    """Reads a VASP CONTCAR file to get structure information."""
    print(f"Reading structure data from {filepath}...")
    with open(filepath, 'r') as f:
        lines = f.readlines()

    comment = lines[0].strip()
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
    
    positions_lines = [lines[i] for i in range(pos_start_idx, pos_start_idx + total_atoms)]
    positions = np.array([list(map(float, line.split()[:3])) for line in positions_lines])

    # Check if coordinates are Cartesian or Direct
    if lines[pos_start_idx - 1].strip().lower().startswith('c'):
        positions_cart = positions
        positions_frac = positions @ np.linalg.inv(lattice_vectors)
    else: # Assume Direct (fractional)
        positions_frac = positions
        positions_cart = positions_frac @ lattice_vectors

    # Build a full list of atom symbols
    atom_symbols_full = []
    for sym, count in zip(atom_symbols, atom_counts):
        atom_symbols_full.extend([sym] * count)

    return {
        "lattice_vectors": lattice_vectors,
        "atom_symbols": atom_symbols,
        "atom_symbols_full": atom_symbols_full,
        "atom_counts": atom_counts,
        "total_atoms": total_atoms,
        "positions_cart": positions_cart,
        "positions_frac": positions_frac,
        "header_lines": header
    }

def read_band_yaml(filepath="band.yaml"):
    """Parses a Phonopy band.yaml file to get phonon modes."""
    print(f"Reading phonon modes from {filepath}...")
    try:
        with open(filepath, 'r') as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: '{filepath}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing YAML file: {e}")
        sys.exit(1)
    
    try:
        n_atoms = data['natom']
        phonons = data['phonon'][0]['band']
        n_modes = len(phonons)
        masses = np.array([atom['mass'] for atom in data['points']])
        
        frequencies = np.array([p['frequency'] for p in phonons])
        eigenvectors_raw = np.array([p['eigenvector'] for p in phonons])
        
        eigenvectors = eigenvectors_raw[..., 0].reshape(n_modes, n_modes)
        
        masses_expanded = np.repeat(masses, 3)
        eigendisplacements = eigenvectors / np.sqrt(masses_expanded)
        eigendisplacements = eigendisplacements.reshape(n_modes, n_atoms, 3)
        
        return frequencies, eigendisplacements, masses, n_atoms
    except KeyError as e:
        print(f"Error: Missing key {e} in band.yaml. The file structure may be different.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred reading band.yaml: {e}")
        sys.exit(1)


def write_vmd_script(filename, positions, displacements, l_cylinder, l_cone):
    """Writes the VMD drawing commands for a single mode."""
    with open(filename, 'w') as f:
        f.write("draw color blue\n")
        for i in range(len(positions)):
            pos, disp = positions[i], displacements[i]
            p1 = pos
            p2 = pos + l_cylinder * disp
            p3 = p2 + l_cone * disp
            
            f.write(f"draw cylinder {{{p1[0]:.6f} {p1[1]:.6f} {p1[2]:.6f}}} {{{p2[0]:.6f} {p2[1]:.6f} {p2[2]:.6f}}} radius 0.1 resolution 30 filled yes\n")
            f.write(f"draw cone {{{p2[0]:.6f} {p2[1]:.6f} {p2[2]:.6f}}} {{{p3[0]:.6f} {p3[1]:.6f} {p3[2]:.6f}}} radius 0.2 resolution 30\n")

def write_vesta_file(filename, template_content, displacements, n_atoms, scale_factor, freq_cm1):
    """
    Writes a new .vesta file by injecting vector data into a template
    based on a working example.
    """
    
    insert_pos = -1
    markers = ['SPLAN', 'VECTR', 'VECTT', 'BOUNDS', 'MODGS', 'END']
    
    for marker in markers:
        pos = template_content.find(marker)
        if pos != -1:
            if insert_pos == -1 or pos < insert_pos:
                insert_pos = pos
    
    if insert_pos == -1:
         print(f"Error: The template file is missing a 'SPLAN' or 'BOUNDS' marker.")
         print("Please open the file in VESTA and re-save it to ensure it has all default sections.")
         return False
    else:
        header = template_content[:insert_pos]
        footer = template_content[insert_pos:]

    header = re.sub(r'VECTR.*', '', header, flags=re.DOTALL)
    header = re.sub(r'VECTT.*', '', header, flags=re.DOTALL)
    header = header.strip()
    
    content = [header]
    content.append("\n# Phonon Mode Vector Data")
    content.append(f"# Frequency: {freq_cm1:.2f} cm-1")
    
    content.append('VECTR')
    scaled_disps = displacements * scale_factor
    for i in range(n_atoms):
        vec = scaled_disps[i]
        content.append(f"  {i+1:3d}  {vec[0]:10.5f}{vec[1]:10.5f}{vec[2]:10.5f} 0")
        content.append(f"   {i+1:3d}   0    0    0    0")
        content.append(" 0 0 0 0 0")
    content.append(" 0 0 0 0 0")
    
    content.append('VECTT')
    for i in range(n_atoms):
        content.append(f"   {i+1:3d} 0.250 255   0   0 0") # Style: 0.25Å radius, Red
    content.append(" 0 0 0 0 0")
    
    content.append(footer)
    
    with open(filename, 'w') as f:
        f.write("\n".join(content))
    return True

def run_visualization():
    """Main function to generate visualization scripts for phonon modes."""
    required_files = ["CONTCAR", "band.yaml"]
    for f in required_files:
        if not os.path.exists(f):
            print(f"***** {f} not found *****")
            sys.exit(1)

    print("Reading CONTCAR and band.yaml...")
    structure = read_contcar()
    frequencies, eigendisps, masses, n_atoms = read_band_yaml()
    n_modes = len(frequencies)

    try:
        factor = float(input("Please input the arrow scaling factor (e.g., 0.5): "))
        output_format = input("Select output format [V]MD, V[E]STA, or [B]oth: ").lower()
    except ValueError:
        print("Invalid factor. Exiting.")
        sys.exit(1)
        
    if output_format not in ['v', 'e', 'b']:
        print("Invalid format choice. Defaulting to Both.")
        output_format = 'b'
        
    write_vmd = output_format in ['v', 'b']
    write_vesta = output_format in ['e', 'b']
    
    vesta_template_content = ""
    if write_vesta:
        template_name = input("Enter the name of your template .vesta file (e.g., MoS2.vesta): ")
        if not os.path.exists(template_name) or not template_name.endswith(".vesta"):
            print(f"Error: Template file '{template_name}' not found or is not a .vesta file.")
            print("Please create one by opening your CONTCAR in VESTA and using 'File -> Save As...'.")
            write_vesta = False
        else:
            with open(template_name, 'r') as f:
                vesta_template_content = f.read()

    # Scale factor
    factor *= np.sqrt(np.max(masses))
    thz_to_cm1 = 33.35641
    
    # VMD arrow scaling
    l_cylinder = 4.0 * factor
    l_cone = 1.5 * factor
    # VESTA arrow scaling (total length)
    vesta_scale = l_cylinder + l_cone

    print("Generating mode visualization files...")
    os.makedirs("VMD_MODES", exist_ok=True)
    os.makedirs("VESTA_MODES", exist_ok=True)
    
    for i in range(n_modes):
        freq_cm1 = frequencies[i] * thz_to_cm1
        
        if write_vmd:
            filename_vmd = os.path.join("VMD_MODES", f"mode_{i+1:03d}.vmd")
            write_vmd_script(filename_vmd, structure["positions_cart"], eigendisps[i], l_cylinder, l_cone)
        
        if write_vesta:
            filename_vesta = os.path.join("VESTA_MODES", f"mode_{i+1:03d}_({freq_cm1:.1f}cm-1).vesta")
            write_vesta_file(
                filename_vesta, 
                vesta_template_content, 
                eigendisps[i], 
                n_atoms, 
                factor,
                freq_cm1 
            )

    print(f"\nGenerated {n_modes} file(s) for the primitive cell.")
    if write_vmd:
        print(f"VMD files are in 'VMD_MODES/'")
    if write_vesta:
        print(f"VESTA files are in 'VESTA_MODES/'")

    if not os.path.exists("all_modes.txt"):
        print("Writing all_modes.txt...")
        with open("all_modes.txt", "w") as f:
            f.write("# mode  freq(cm-1)\n")
            for i in range(n_modes):
                f.write(f"{i+1:4d}   {frequencies[i] * thz_to_cm1:10.4f}\n")

    supercell_choice = input("Visualize in a supercell? (yes/no): ").lower()
    if supercell_choice == 'yes':
        try:
            supercell_dim_str = input("Please input the dimension of supercell (e.g., 2 2 1): ")
            supercell_dim = np.array(supercell_dim_str.split(), dtype=int)
        except (ValueError, IndexError):
            print("Invalid supercell dimension. Exiting.")
            sys.exit(1)

        print("Generating supercell and visualization files...")
        
        supercell_positions = []
        prim_pos = structure["positions_cart"]
        lv = structure["lattice_vectors"]
        
        for i in range(supercell_dim[0]):
            for j in range(supercell_dim[1]):
                for k in range(supercell_dim[2]):
                    translation = i * lv[0] + j * lv[1] + k * lv[2]
                    supercell_positions.append(prim_pos + translation)
        
        supercell_positions = np.concatenate(supercell_positions, axis=0)
        
        supercell_lvs = structure["lattice_vectors"] * supercell_dim[:, np.newaxis]
        
        supercell_atom_symbols_full = []
        for i in range(np.prod(supercell_dim)):
            supercell_atom_symbols_full.extend(structure["atom_symbols_full"])

        structure_supercell = {
            "lattice_vectors": supercell_lvs,
            "atom_symbols": structure["atom_symbols"],
            "atom_symbols_full": supercell_atom_symbols_full,
            "atom_counts": structure["atom_counts"] * np.prod(supercell_dim),
            "total_atoms": structure["total_atoms"] * np.prod(supercell_dim),
            "positions_cart": supercell_positions,
            "positions_frac": supercell_positions @ np.linalg.inv(supercell_lvs)
        }
        
        with open("POSCAR_supercell", "w") as f:
            f.write("Supercell for visualization\n")
            f.write("1.0\n")
            for vec in structure_supercell["lattice_vectors"]:
                f.write(f"   {vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}\n")
            supercell_atom_counts = structure["atom_counts"] * np.prod(supercell_dim)
            f.write("  " + " ".join(structure_supercell["atom_symbols"]) + "\n")
            f.write("  " + " ".join(map(str, supercell_atom_counts)) + "\n")
            f.write("Cartesian\n")
            for pos in structure_supercell["positions_cart"]:
                f.write(f"   {pos[0]:.9f} {pos[1]:.9f} {pos[2]:.9f}\n")
                
        for i in range(n_modes):
            supercell_disps = np.tile(eigendisps[i], (np.prod(supercell_dim), 1))
            freq_cm1 = frequencies[i] * thz_to_cm1
            
            if write_vmd:
                filename_vmd = os.path.join("VMD_MODES", f"mode_super_{i+1:03d}.vmd")
                write_vmd_script(filename_vmd, structure_supercell["positions_cart"], supercell_disps, l_cylinder, l_cone)
            if write_vesta:
                filename_vesta = os.path.join("VESTA_MODES", f"mode_super_{i+1:03d}_({freq_cm1:.1f}cm-1).vesta")
                # We need to re-read the template content because the original template is for the primitive cell
                write_vesta_file(
                    filename_vesta,
                    vesta_template_content, # Note: This is still the primitive cell template
                    supercell_disps,
                    structure_supercell["total_atoms"],
                    factor,
                    freq_cm1
                )
        
        print(f"Generated {n_modes} file(s) for the supercell.")
            
    print("\nvisualize_modes.py finished successfully.")

if __name__ == "__main__":
    run_visualization()