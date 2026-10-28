from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import yaml
import json
import hashlib

from spectropy_structure import read_structure


@dataclass(frozen=True)
class GammaModes:
    frequencies_thz: np.ndarray
    eigendisplacements: np.ndarray
    masses: np.ndarray
    natoms: int


def read_gamma_modes(path: str | Path = "band.yaml", complex_eigenvectors: bool = False) -> GammaModes:
    with Path(path).open() as stream:
        data = yaml.load(stream, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    natoms = data["natom"]
    if complex_eigenvectors and not np.allclose(data["phonon"][0]["q-position"], 0, atol=1e-8):
        raise ValueError("Complex mode diagnostics require Gamma-point modes as the first band.yaml entry")
    bands = data["phonon"][0]["band"]
    frequencies = np.array([band["frequency"] for band in bands])
    pairs = np.array([band["eigenvector"] for band in bands])
    eigenvectors = pairs[..., 0] + 1j * pairs[..., 1] if complex_eigenvectors else pairs[..., 0]
    masses = np.array([point["mass"] for point in data["points"]])
    nmodes = natoms * 3
    eigendisplacements = eigenvectors.reshape(nmodes, nmodes) / np.sqrt(np.repeat(masses, 3))
    return GammaModes(frequencies, eigendisplacements.reshape(nmodes, natoms, 3), masses, natoms)


def align_mode_atoms(phonons, structure, band_path="band.yaml", tolerance=1e-3):
    """Match band.yaml atoms to derivative atoms, allowing a common origin shift.

    No structure or phonon file is changed. Lattices must already have the
    same Cartesian orientation. Return aligned modes and mapping provenance.
    """
    lines = []
    in_phonons = False
    with Path(band_path).open() as stream:
        for line in stream:
            if line.startswith("phonon:"):
                in_phonons = True
                continue
            if in_phonons and line and not line[0].isspace() and not line.startswith(("-", "#")) and ":" in line:
                in_phonons = False
            if not in_phonons:
                lines.append(line)
    metadata = yaml.safe_load("".join(lines))
    if not np.allclose(metadata["lattice"], structure.lattice, atol=1e-6, rtol=0):
        raise ValueError("Phonon and derivative lattices differ; explicit Cartesian frame alignment required")
    positions = np.asarray([p["coordinates"] for p in metadata["points"]])
    symbols = np.asarray([p["symbol"] for p in metadata["points"]])
    reference_symbols = np.asarray(structure.atom_symbols)
    if len(positions) != structure.natoms or phonons.natoms != structure.natoms:
        raise ValueError("Phonon and derivative atom counts differ")
    if not np.isfinite(positions).all() or not np.isfinite(phonons.masses).all() or np.any(phonons.masses <= 0):
        raise ValueError("Invalid phonon coordinates or masses")
    def match(shift):
        delta = positions[:, None, :] + shift - structure.fractional_positions[None, :, :]
        delta -= np.rint(delta)
        distances = np.linalg.norm(delta @ structure.lattice, axis=2)
        distances[symbols[:, None] != reference_symbols[None, :]] = np.inf
        mapping = distances.argmin(axis=1)
        maximum = distances[np.arange(len(mapping)), mapping].max()
        return maximum, mapping
    maximum, mapping = match(np.zeros(3))
    shift = np.zeros(3)
    if maximum > tolerance or len(set(mapping)) != len(mapping):
        candidates = []
        for target in np.where(reference_symbols == symbols[0])[0]:
            trial_shift = structure.fractional_positions[target] - positions[0]
            trial_max, trial_mapping = match(trial_shift)
            if trial_max <= tolerance and len(set(trial_mapping)) == len(trial_mapping):
                candidates.append((trial_max, trial_shift, trial_mapping))
        if len(candidates) != 1:
            raise ValueError("No unique phonon-to-derivative atom alignment; provide consistently ordered input files")
        maximum, shift, mapping = candidates[0]
    eig = np.empty_like(phonons.eigendisplacements)
    eig[:, mapping] = phonons.eigendisplacements
    masses = np.empty_like(phonons.masses)
    masses[mapping] = phonons.masses
    aligned = GammaModes(phonons.frequencies_thz, eig, masses, phonons.natoms)
    return aligned, mapping, shift, maximum


def read_registered_modes(path="band.yaml", reference=None, provenance=None):
    """Require geometry registration before any atomic-response contraction."""
    path = Path(path)
    if reference is None:
        reference = path.parent / "ref_poscar.vasp"
        if not reference.is_file():
            reference = path.parent / "CONTCAR"
    reference = Path(reference)
    if not reference.is_file():
        raise FileNotFoundError("Atom registration requires ref_poscar.vasp or an explicit derivative reference")
    phonons = read_gamma_modes(path, complex_eigenvectors=True)
    aligned, mapping, shift, error = align_mode_atoms(phonons, read_structure(reference), path)
    if provenance is not None:
        report = dict(band=str(path.resolve()), reference=str(reference.resolve()),
            band_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
            phonon_to_derivative_1based=(mapping + 1).tolist(),
            fractional_origin_shift=shift.tolist(), maximum_error_angstrom=float(error),
            tolerance_angstrom=1e-3)
        Path(provenance).write_text(json.dumps(report, indent=2) + "\n")
    return aligned


def read_irreps(path: str | Path = "irreps.yaml") -> list[str] | None:
    file_path = Path(path)
    if not file_path.is_file():
        return None
    labels: dict[int, str] = {}
    for mode in yaml.safe_load(file_path.read_text())["normal_modes"]:
        for band_index in mode["band_indices"]:
            labels[band_index] = mode["ir_label"]
    return [labels.get(index + 1, "---") for index in range(len(labels))]
