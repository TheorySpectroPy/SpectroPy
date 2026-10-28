import os
import shutil
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

def read_diel_from_outcar(path, target_frequency):
    """Read density-density optics at an exact photon energy, retaining digits."""
    with open(path) as stream:
        lines = stream.read().splitlines()
    tensors = []
    for kind in ('REAL', 'IMAGINARY'):
        start = next((i for i, line in enumerate(lines)
                      if 'DIELECTRIC FUNCTION' in line and kind in line
                      and 'density-density' in line), None)
        if start is None:
            raise ValueError(f'Missing {kind} density-density dielectric table in {path}')
        rows = []
        for line in lines[start + 3:]:
            parts = line.split()
            if len(parts) != 7:
                break
            try:
                rows.append([float(value) for value in parts])
            except ValueError:
                break
        values = np.asarray(rows)
        if values.ndim != 2 or not np.isfinite(values).all():
            raise ValueError(f'Invalid {kind} dielectric table in {path}')
        row = values[np.argmin(abs(values[:, 0] - target_frequency))]
        if abs(row[0] - target_frequency) > 1e-8:
            raise ValueError(f'{target_frequency:g} eV is not present on the OUTCAR grid in {path}; '
                             'choose an available energy or dielectric_source: vasprun.xml')
        x, y, z, xy, yz, zx = row[1:]
        tensors.append(np.array([[x, xy, zx], [xy, y, yz], [zx, yz, z]]))
    return tuple(tensors)


def read_dielectric(path, target_frequency, source='outcar'):
    if source == 'outcar':
        return read_diel_from_outcar(path, target_frequency)
    if source == 'vasprun.xml':
        return read_diel_from_xml(path, target_frequency)
    raise ValueError('dielectric_source must be outcar or vasprun.xml')

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


def collect_vasprun_files(displacements, destination=None, source="vasprun.xml"):
    from spectropy_paths import displacement_root
    root = displacement_root()
    destination = root / 'vasprun' if destination is None else destination
    os.makedirs(destination, exist_ok=True)
    for displacement in displacements:
        filename = os.path.join(root, f"ra_pos_{displacement.name}", 'OUTCAR' if source == 'outcar' else 'vasprun.xml')
        try:
            suffix = 'OUTCAR' if source == 'outcar' else 'xml'
            shutil.copyfile(filename, os.path.join(destination, f"{displacement.name}.{suffix}"))
        except FileNotFoundError:
            print(f"Error: Source file not found: {filename}")
            print("Please ensure all VASP calculations have finished.")
            sys.exit(1)


def read_displaced_dielectrics(displacements, lattice, target_frequency, directory=None, source="vasprun.xml"):
    from spectropy_paths import displacement_root
    if directory is None:
        directory = displacement_root() / 'vasprun'
    measured = defaultdict(list)
    for displacement in displacements:
        suffix = 'OUTCAR' if source == 'outcar' else 'xml'
        filename = os.path.join(directory, f"{displacement.name}.{suffix}")
        real, imaginary = read_dielectric(filename, target_frequency, source)
        if real is None:
            sys.exit(1)
        direction = displacement.fractional_vector @ lattice
        measured[displacement.atom_index].append((direction, real, imaginary))
    return dict(measured)
