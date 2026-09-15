import os
import shutil
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

def read_diel_from_xml(xml_path, target_frequency):
    try:
        root = ET.parse(xml_path).getroot()
    except (ET.ParseError, FileNotFoundError):
        print(f"Error: Could not parse or find {xml_path}")
        return None, None
    try:
        dielectric = root.find("./calculation/dielectricfunction")
        imaginary = dielectric.find("./imag/array/set")
        real = dielectric.find("./real/array/set")
        frequencies = np.array([float(row.text.split()[0]) for row in imaginary.findall("r")])
        real_values = np.array([list(map(float, row.text.split()[1:])) for row in real.findall("r")])
        imaginary_values = np.array([list(map(float, row.text.split()[1:])) for row in imaginary.findall("r")])
        index = np.argmin(np.abs(frequencies - target_frequency))
        real_components = real_values[index]
        imaginary_components = imaginary_values[index]
        real_tensor = np.array([[real_components[0], real_components[3], real_components[5]], [real_components[3], real_components[1], real_components[4]], [real_components[5], real_components[4], real_components[2]]])
        imaginary_tensor = np.array([[imaginary_components[0], imaginary_components[3], imaginary_components[5]], [imaginary_components[3], imaginary_components[1], imaginary_components[4]], [imaginary_components[5], imaginary_components[4], imaginary_components[2]]])
        return real_tensor, imaginary_tensor
    except AttributeError:
        print(f"Error: Could not find dielectric data in {xml_path}")
        return None, None


def collect_vasprun_files(displacements, destination="vasprun"):
    os.makedirs(destination, exist_ok=True)
    for displacement in displacements:
        source = os.path.join(f"ra_pos_{displacement.name}", "vasprun.xml")
        try:
            shutil.copyfile(source, os.path.join(destination, f"{displacement.name}.xml"))
        except FileNotFoundError:
            print(f"Error: Source file not found: {source}")
            print("Please ensure all VASP calculations have finished.")
            sys.exit(1)


def read_displaced_dielectrics(displacements, lattice, target_frequency, directory="vasprun"):
    measured = defaultdict(list)
    for displacement in displacements:
        filename = os.path.join(directory, f"{displacement.name}.xml")
        real, imaginary = read_diel_from_xml(filename, target_frequency)
        if real is None:
            sys.exit(1)
        direction = displacement.fractional_vector @ lattice
        measured[displacement.atom_index].append((direction, real, imaginary))
    return dict(measured)
