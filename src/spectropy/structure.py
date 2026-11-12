from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class Structure:
    comment: str
    lattice: np.ndarray
    symbols: tuple[str, ...]
    counts: np.ndarray
    fractional_positions: np.ndarray

    @property
    def natoms(self) -> int:
        return int(self.counts.sum())

    @property
    def cartesian_positions(self) -> np.ndarray:
        return self.fractional_positions @ self.lattice

    @property
    def atom_symbols(self) -> list[str]:
        return [symbol for symbol, count in zip(self.symbols, self.counts) for _ in range(int(count))]


def read_structure(path: str | Path = "CONTCAR") -> Structure:
    lines = Path(path).read_text().splitlines()
    if len(lines) < 8:
        raise ValueError(f"Incomplete VASP structure: {path}")
    scale = float(lines[1])
    lattice = np.array([list(map(float, line.split()[:3])) for line in lines[2:5]]) * scale
    try:
        counts = np.array(lines[6].split(), dtype=int)
        symbols = tuple(lines[5].split())
        coordinate_line = 7
    except ValueError:
        counts = np.array(lines[5].split(), dtype=int)
        symbols = tuple(f"Elem{index + 1}" for index in range(len(counts)))
        coordinate_line = 6
    if lines[coordinate_line].lower().startswith("s"):
        coordinate_line += 1
    coordinates = np.array([
        list(map(float, line.split()[:3]))
        for line in lines[coordinate_line + 1: coordinate_line + 1 + int(counts.sum())]
    ])
    if lines[coordinate_line].lower().startswith(("c", "k")):
        coordinates = coordinates @ np.linalg.inv(lattice)
    return Structure(lines[0], lattice, symbols, counts, coordinates)


def write_structure(path: str | Path, structure: Structure, fractional_positions: np.ndarray | None = None) -> None:
    positions = structure.fractional_positions if fractional_positions is None else fractional_positions
    with open(path, "w") as output:
        output.write(f"{structure.comment}\n1.0\n")
        for vector in structure.lattice:
            output.write(f"{vector[0]:21.16f} {vector[1]:21.16f} {vector[2]:21.16f}\n")
        output.write(" ".join(structure.symbols) + "\n")
        output.write(" ".join(map(str, structure.counts)) + "\nDirect\n")
        for position in positions:
            output.write(f"{position[0]:19.16f} {position[1]:19.16f} {position[2]:19.16f}\n")
