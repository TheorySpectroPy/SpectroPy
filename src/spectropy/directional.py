"""Measured-direction raw output, with explicit physical units."""
from pathlib import Path
import numpy as np
from .reconstruction import fractional_finite_differences


def directional_document(energy, structure, measurements, sites,
                         scheme='central', equilibrium=None):
    lattice = np.asarray(structure.lattice)
    factor = abs(np.linalg.det(lattice)) / (4*np.pi)
    samples = []
    eq = None if equilibrium is None else equilibrium[0] + 1j*equilibrium[1]
    for atom, values in measurements.items():
        rotations = sites.get(atom, {'rotations': np.eye(3)[None]})['rotations']
        measured = [(v, real + 1j*imag) for v, real, imag in values]
        for index, (unit, response, length, generated) in enumerate(
                fractional_finite_differences(measured, rotations, lattice,
                                             difference_scheme=scheme, equilibrium=eq), 1):
            cart = unit @ lattice
            scale = np.linalg.norm(cart)
            response = response * factor / scale
            samples.append(dict(atom=int(atom), direction_index=index,
                unit_cartesian=(cart/scale).tolist(), unit_fractional=unit.tolist(),
                displacement_length_angstrom=float(length*scale),
                negative_from_symmetry=generated,
                calculated_signs=1 if generated or scheme == 'forward' else 2,
                real=response.real.tolist(), imaginary=response.imag.tolist()))
    if not samples:
        raise ValueError('No measured directional responses')
    return dict(energy_eV=float(energy),
        response_units='angstrom^2', dielectric_basis='cartesian',
        displacement_basis='measured unit directions',
        volume_factor_angstrom3=float(factor), difference_scheme=scheme,
        lattice_angstrom=lattice.tolist(), natoms=int(structure.natoms),
        samples=samples)


def save_directional(path, document):
    """Write measured-direction 3x3 matrices, not inferred Cartesian slices."""
    def numbers(vector):
        return " ".join(f"{value:.16e}" for value in vector)

    lattice = np.asarray(document["lattice_angstrom"], dtype=float)
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
        raise ValueError("Invalid directional derivative lattice")
    if document["response_units"] != "angstrom^2":
        raise ValueError("Unexpected directional derivative units")
    if not np.isfinite([document["energy_eV"], document["volume_factor_angstrom3"]]).all():
        raise ValueError("Nonfinite directional derivative metadata")
    lines = [
        f"! Raw directional dielectric derivatives at {document['energy_eV']:.16e} eV",
        "! Matrix basis: Cartesian x,y,z; displacement basis: measured unit directions.",
        "! Only measured directions are listed; no denoising or equivalent-atom reconstruction.",
        "! Units: angstrom^2; responses include Vcell/(4*pi).",
        f"! Volume factor (angstrom^3): {document['volume_factor_angstrom3']:.16e}",
        f"! Difference scheme: {document['difference_scheme']}",
        f"! Total atoms: {document['natoms']}",
        "! Unit cell matrix (row vectors, angstrom):",
        *["!   " + numbers(row) for row in lattice],
    ]
    seen = set()
    for sample in document["samples"]:
        key = (sample["atom"], sample["direction_index"])
        if key in seen or not 1 <= key[0] <= document["natoms"] or key[1] < 1:
            raise ValueError("Duplicate or invalid directional sample")
        seen.add(key)
        arrays = {}
        for name, shape in (("unit_cartesian", (3,)), ("unit_fractional", (3,)),
                            ("real", (3, 3)), ("imaginary", (3, 3))):
            arrays[name] = np.asarray(sample[name], dtype=float)
            if arrays[name].shape != shape or not np.isfinite(arrays[name]).all():
                raise ValueError(f"Invalid directional derivative {name}")
        step = sample["displacement_length_angstrom"]
        if not np.isfinite(step) or step <= 0:
            raise ValueError("Displacement length must be finite and positive")
        lines.extend([
            "",
            f"Atom {key[0]} Direction {key[1]}",
            "! Cartesian unit direction: " + numbers(arrays["unit_cartesian"]),
            "! Fractional direction: " + numbers(arrays["unit_fractional"]),
            f"! Displacement length (angstrom): {step:.16e}",
            f"! Calculated signs: {sample['calculated_signs']}",
            "! Negative from symmetry: " + str(bool(sample["negative_from_symmetry"])).lower(),
            "! Real Part",
            *[numbers(row) for row in arrays["real"]],
            "! Imaginary Part",
            *[numbers(row) for row in arrays["imaginary"]],
        ])
    if not seen:
        raise ValueError("No measured directional responses")
    Path(path).write_text("\n".join(lines) + "\n")
