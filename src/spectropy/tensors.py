"""Final mode tensors: in-memory atomic contraction and lossless text I/O."""
from pathlib import Path
import re

import numpy as np

from .config import DEFAULT_INPUT, resolve_settings
from .phonons import read_registered_modes, read_irreps
from .structure import read_structure


THZ_TO_CM1 = 33.35641


def write_raman_tensors(path, frequencies, tensors, labels=None, mode_indices=None):
    frequencies = np.asarray(frequencies, dtype=float)
    tensors = np.asarray(tensors, dtype=complex)
    count = len(frequencies)
    indices = np.arange(1, count + 1) if mode_indices is None else np.asarray(mode_indices)
    if (frequencies.shape != (count,) or tensors.shape != (count, 3, 3)
            or not np.isfinite(frequencies).all() or not np.isfinite(tensors).all()
            or indices.shape != (count,) or np.any(indices < 1)
            or not np.issubdtype(indices.dtype, np.integer) or len(set(indices)) != count):
        raise ValueError("Invalid mode tensor data")
    labels = [labels[i] if labels is not None and i < len(labels) else "---" for i in range(count)]
    if any(not isinstance(label, str) or len(label.split()) != 1 for label in labels):
        raise ValueError("Mode labels must be single nonempty tokens")
    with Path(path).open("w") as stream:
        stream.write("# Mode Freq(THz) Freq(cm-1) Irrep. Raman Tensor (Real + i*Imaginary)\n")
        for index, frequency, tensor, label in zip(indices, frequencies, tensors, labels):
            stream.write(f"{index} {frequency:.16e} {frequency * THZ_TO_CM1:.16e} {label}\n")
            for row in tensor:
                stream.write(" ".join(f"{z.real:.16e}{z.imag:+.16e}j" for z in row) + "\n")
            stream.write("\n")


def read_raman_tensors(path):
    """Return mode indices, frequencies, tensors, and labels in file order.

    Mode indices remain 1-based and may describe only a subset of modes.
    No transpose symmetry or rounding is imposed.
    """
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    number = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    entry = rf"([+-]?{number})\s*([+-]{number})j"
    row_pattern = re.compile(rf"\s*{entry}\s+{entry}\s+{entry}\s*")
    indices, frequencies, tensors, labels = [], [], [], []
    try:
        if not lines or len(lines) % 4:
            raise ValueError("incomplete tensor blocks")
        for offset in range(0, len(lines), 4):
            header = lines[offset].split()
            if len(header) not in (3, 4):
                raise ValueError("invalid mode header")
            index, frequency, frequency_cm1 = int(header[0]), float(header[1]), float(header[2])
            if index < 1 or index in indices or not np.isfinite([frequency, frequency_cm1]).all():
                raise ValueError("invalid or duplicate mode")
            tensor = []
            for line in lines[offset + 1:offset + 4]:
                match = row_pattern.fullmatch(line)
                if match is None:
                    raise ValueError("expected three complex entries per row")
                values = match.groups()
                tensor.append([complex(float(values[i]), float(values[i + 1])) for i in (0, 2, 4)])
            indices.append(index)
            frequencies.append(frequency)
            tensors.append(tensor)
            labels.append(header[3] if len(header) == 4 else "---")
        if not np.isfinite(tensors).all():
            raise ValueError("nonfinite tensor")
    except (ValueError, OverflowError) as error:
        raise ValueError(f"Invalid Raman tensor file {path}: {error}") from error
    return np.asarray(indices), np.asarray(frequencies), np.asarray(tensors), labels


def contract_atomic_derivatives(structure, field, phonons):
    """Contract complex D[atom,alpha,beta,i] in 1/Angstrom with aligned modes."""
    if phonons.natoms != structure.natoms:
        raise ValueError("Incomplete atomic field or inconsistent phonon atom count")
    field = np.asarray(field, dtype=complex)
    if field.shape != (structure.natoms, 3, 3, 3) or not np.isfinite(field).all():
        raise ValueError("Invalid atomic derivative field")
    # Eigen-displacements already include inverse square-root atomic masses.
    return abs(np.linalg.det(structure.lattice)) / (4 * np.pi) * np.einsum(
        "jabi,mji->mab", field, phonons.eigendisplacements)


def run_tensors(input_path=DEFAULT_INPUT, use_saved_basis=None,
                equilibrium_vasprun=None, *, settings=None):
    from .derivatives import iter_atomic_derivatives, write_atomic_outputs

    # Validate options and atom registration before processing large XML files.
    settings = resolve_settings(settings, input_path, saved_basis=use_saved_basis,
                                equilibrium_file=equilibrium_vasprun)
    structure = read_structure('CONTCAR')
    phonons = read_registered_modes(reference="CONTCAR")
    labels = read_irreps()

    for energy, field, directional, measurements, equilibrium in iter_atomic_derivatives(settings, structure):
        tensors = contract_atomic_derivatives(structure, field, phonons)
        path = f"Raman_tensor_{energy:.2f}eV"
        write_raman_tensors(path, phonons.frequencies_thz, tensors, labels)
        print(f"Wrote {path}")
        if settings.per_atom_output:
            write_atomic_outputs(energy, structure, field, directional, measurements, equilibrium, settings)
