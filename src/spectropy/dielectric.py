import os
import shutil
import sys
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

def _xml_only(source):
    if source != 'vasprun.xml':
        raise ValueError('dielectric_source must be vasprun.xml; other VASP dielectric readers were removed')


def read_dielectric(path, target_frequency, source='vasprun.xml'):
    _xml_only(source)
    return read_diel_from_xml(path, target_frequency)

def read_diel_from_xml(xml_path, target_frequency):
    return read_dielectrics_xml(xml_path, [target_frequency])[0]


def read_dielectrics_xml(xml_path, energies):
    """Read density-density XML once, preserving every printed component.

    Sample the nearest tabulated photon energy, without interpolation. Missing
    off-diagonal transpose entries are filled by reciprocity with a warning.
    """
    root = ET.parse(xml_path).getroot()
    blocks = root.findall('.//dielectricfunction')
    selected = [b for b in blocks if b.get('comment', '').strip() == 'density-density']
    if not selected:
        selected = [b for b in blocks if not b.get('comment', '').strip()]
    if len(selected) != 1:
        raise ValueError(f'Expected one density-density dielectric block in {xml_path}')
    grids=[];tensors=[];filled=[]
    for kind in ('real', 'imag'):
        array = selected[0].find(f'{kind}/array')
        if array is None:
            raise ValueError(f'Missing {kind} dielectric data in {xml_path}')
        fields = [(f.text or '').strip().lower() for f in array.findall('field')]
        expected = {'energy'} | {a+b for a in 'xyz' for b in 'xyz'}
        if (len(fields) != len(set(fields)) or not set(fields) <= expected or
                not {'energy', 'xx', 'yy', 'zz'} <= set(fields) or
                any(a+b not in fields and b+a not in fields for a,b in [('x','y'),('x','z'),('y','z')])):
            raise ValueError(f'{xml_path}: {kind} requires unique named fields: energy, xx, yy, zz, '
                             'and at least one component of each xy/yx, xz/zx, yz/zy pair')
        rows = np.asarray([[float(v) for v in (r.text or '').split()] for r in array.findall('.//r')])
        if rows.ndim != 2 or rows.shape[1] != len(fields) or not len(rows) or not np.isfinite(rows).all():
            raise ValueError(f'Invalid {kind} dielectric rows in {xml_path}')
        grid = rows[:, fields.index('energy')]
        if np.any(np.diff(grid) <= 0):
            raise ValueError(f'Non-increasing dielectric grid in {xml_path}')
        grids.append(grid)
        columns = []
        for a in 'xyz':
            for b in 'xyz':
                component = a+b
                if component not in fields:
                    filled.append(f'{kind} {component}={b+a}')
                    component = b+a
                columns.append(fields.index(component))
        tensors.append(rows[:, columns].reshape(-1,3,3))
    if not np.array_equal(grids[0], grids[1]):
        raise ValueError(f'Real and imaginary dielectric grids differ in {xml_path}')
    if filled:
        warnings.warn(f'{xml_path}: not all nine dielectric components are printed. '
                      f'Assuming reciprocity to fill {", ".join(filled)} (no complex conjugation). '
                      'These entries are inferred, not independently calculated; '
                      'nonreciprocal response cannot be recovered from missing components.',
                      RuntimeWarning, stacklevel=2)
    result=[]
    for energy in energies:
        if not np.isfinite(energy) or energy < grids[0][0] or energy > grids[0][-1]:
            raise ValueError(f'Photon energy {energy} outside dielectric grid in {xml_path}')
        i = np.argmin(np.abs(grids[0] - energy))
        result.append((tensors[0][i],tensors[1][i]))
    return result


def collect_vasprun_files(displacements, destination=None, source="vasprun.xml"):
    _xml_only(source)
    from .paths import displacement_root
    root = displacement_root()
    destination = root / 'vasprun' if destination is None else destination
    os.makedirs(destination, exist_ok=True)
    for displacement in displacements:
        filename = os.path.join(root, f"ra_pos_{displacement.name}", 'vasprun.xml')
        try:
            shutil.copyfile(filename, os.path.join(destination, f"{displacement.name}.xml"))
        except FileNotFoundError:
            print(f"Error: Source file not found: {filename}")
            print("Please ensure all VASP calculations have finished.")
            sys.exit(1)


def read_displaced_dielectrics(displacements, lattice, target_frequency, directory=None, source="vasprun.xml"):
    _xml_only(source)
    from .paths import displacement_root
    if directory is None:
        directory = displacement_root() / 'vasprun'
    measured = defaultdict(list)
    for displacement in displacements:
        filename = os.path.join(directory, f"{displacement.name}.xml")
        real, imaginary = read_dielectric(filename, target_frequency, source)
        direction = displacement.fractional_vector @ lattice
        measured[displacement.atom_index].append((direction, real, imaginary))
    return dict(measured)


def read_displaced_dielectrics_many(displacements, lattice, energies, directory=None,
                                   source="vasprun.xml"):
    """Read each collected XML once for all energies, preserving request order.

    Keep only detached 3x3 samples: the XML reader's slices otherwise retain
    the full frequency-grid arrays for every displacement in memory.
    No persistent cache or interpolation is introduced. The XML reader warns
    whenever reciprocity is needed to complete missing transpose components.
    """
    _xml_only(source)
    from .paths import displacement_root
    if directory is None:
        directory = displacement_root() / 'vasprun'
    energies = list(energies)
    measured = [defaultdict(list) for _ in energies]
    if not energies:
        return []
    for displacement in displacements:
        filename = os.path.join(directory, f"{displacement.name}.xml")
        samples = read_dielectrics_xml(filename, energies)
        direction = displacement.fractional_vector @ lattice
        for measurements, (real, imaginary) in zip(measured, samples):
            measurements[displacement.atom_index].append(
                (direction.copy(), real.copy(), imaginary.copy()))
    return [dict(measurements) for measurements in measured]
