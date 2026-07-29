from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import yaml


@dataclass(frozen=True)
class GammaModes:
    frequencies_thz: np.ndarray
    eigendisplacements: np.ndarray
    masses: np.ndarray
    natoms: int


def read_gamma_modes(path: str | Path = "band.yaml") -> GammaModes:
    data = yaml.safe_load(Path(path).read_text())
    natoms = data["natom"]
    bands = data["phonon"][0]["band"]
    frequencies = np.array([band["frequency"] for band in bands])
    eigenvectors = np.array([band["eigenvector"] for band in bands])[..., 0]
    masses = np.array([point["mass"] for point in data["points"]])
    nmodes = natoms * 3
    eigendisplacements = eigenvectors.reshape(nmodes, nmodes) / np.sqrt(np.repeat(masses, 3))
    return GammaModes(frequencies, eigendisplacements.reshape(nmodes, natoms, 3), masses, natoms)


def read_irreps(path: str | Path = "irreps.yaml") -> list[str] | None:
    file_path = Path(path)
    if not file_path.is_file():
        return None
    labels: dict[int, str] = {}
    for mode in yaml.safe_load(file_path.read_text())["normal_modes"]:
        for band_index in mode["band_indices"]:
            labels[band_index] = mode["ir_label"]
    return [labels.get(index + 1, "---") for index in range(len(labels))]
