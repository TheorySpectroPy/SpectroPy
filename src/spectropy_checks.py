"""Symmetry consistency diagnostics and polarization-angle scans for Raman tensors."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import yaml

from calculate_spectrum import read_dielectric_derivatives
from spectropy_config import read_settings
from spectropy_phonons import GammaModes, read_gamma_modes, read_irreps, align_mode_atoms, read_registered_modes
from spectropy_structure import read_structure


CM_PER_THZ = 33.35641


def symmetry_operations(structure, symmetry_path="symmetry", tolerance=1e-3):
    """Return atom permutations and Cartesian rotations, preserving species."""
    with Path(symmetry_path).open() as stream:
        data = yaml.load(stream, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    positions = structure.fractional_positions
    species = np.asarray(structure.atom_symbols)
    operations = []
    for entry in data["space_group_operations"]:
        fractional = np.asarray(entry["rotation"], dtype=float)
        rotation = structure.lattice.T @ fractional @ np.linalg.inv(structure.lattice.T)
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6, rtol=0):
            raise ValueError("Symmetry rotation is not orthogonal in Cartesian coordinates")
        moved = positions @ fractional.T + np.asarray(entry["translation"])
        permutation = []
        for index, position in enumerate(moved):
            delta = positions - position
            delta -= np.rint(delta)
            distances = np.linalg.norm(delta @ structure.lattice, axis=1)
            distances[species != species[index]] = np.inf
            target = int(np.argmin(distances))
            if distances[target] > tolerance:
                raise ValueError(f"Symmetry does not map atom {index + 1} onto a matching species")
            permutation.append(target)
        if len(set(permutation)) != structure.natoms:
            raise ValueError("Symmetry atom mapping is not a permutation")
        operations.append((np.asarray(permutation), rotation))
    if not operations:
        raise ValueError("symmetry contains no space_group_operations")
    return operations


def project_derivatives(derivatives, operations):
    """Group-average D[atom, displacement, dielectric_row, dielectric_column]."""
    values = np.asarray(derivatives)
    if values.ndim != 4 or values.shape[1:] != (3, 3, 3):
        raise ValueError("Derivatives must have shape (natoms, 3, 3, 3)")
    if not operations or not np.isfinite(values).all():
        raise ValueError("Finite derivatives and at least one symmetry operation are required")
    result = np.zeros_like(values, dtype=np.result_type(values, float))
    for permutation, rotation in operations:
        result[permutation] += np.einsum(
            "ba,ip,jq,napq->nbij", rotation, rotation, rotation, values, optimize=True
        )
    return result / len(operations)


def mode_tensors(derivatives, eigendisplacements):
    """Contract atomic responses with mass-weighted Gamma eigendisplacements."""
    return np.einsum("naij,mna->mij", derivatives, eigendisplacements, optimize=True)


def projection_residuals(original, projected, epsilon_fraction=1e-12):
    """Return norm residuals with max(original norm, projected norm) regularization."""
    original, projected = np.asarray(original), np.asarray(projected)
    if original.shape != projected.shape or original.ndim != 3 or original.shape[1:] != (3, 3):
        raise ValueError("Mode tensor arrays must share shape (nmodes, 3, 3)")
    if not np.isfinite(epsilon_fraction) or epsilon_fraction <= 0:
        raise ValueError("epsilon_fraction must be finite and positive")
    if not np.isfinite(original).all() or not np.isfinite(projected).all():
        raise ValueError("Mode tensors must be finite")
    a = np.linalg.norm(original, axis=(1, 2))
    b = np.linalg.norm(projected, axis=(1, 2))
    absolute = np.linalg.norm(original - projected, axis=(1, 2))
    epsilon = epsilon_fraction * max(a.max(), b.max(), 1.0)
    return {
        "absolute_residual": absolute,
        "squared_residual": absolute ** 2,
        "relative_residual": absolute / (np.maximum(a, b) + epsilon),
        "original_norm": a,
        "projected_norm": b,
        "epsilon": epsilon,
    }


def mode_representations(phonons, operations, groups, tolerance=1e-5, strict=True):
    """Gamma_pk = <e_p|T(g)e_k> in the supplied mass-weighted mode basis.

    Groups must partition the modes into invariant partner subspaces. No
    textbook partner basis or scalar character is substituted for a multiplet.
    """
    eigenvectors = phonons.eigendisplacements * np.sqrt(phonons.masses)[None, :, None]
    flat = eigenvectors.reshape(len(eigenvectors), -1)
    if not operations or not np.isfinite(flat).all() or not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Finite eigenvectors, operations, and positive tolerance required")
    if not np.allclose(flat.conj() @ flat.T, np.eye(len(flat)), atol=1e-5, rtol=0):
        raise ValueError("Phonon eigenvectors are not orthonormal")
    indices = [index for group in groups for index in group]
    if sorted(indices) != list(range(len(flat))):
        raise ValueError("Irrep partner groups must partition all modes exactly once")
    blocks, leakages = [], []
    for group in groups:
        vectors = eigenvectors[group]
        basis = vectors.reshape(len(group), -1)
        matrices, leakage = [], 0.0
        for permutation, rotation in operations:
            if not np.allclose(phonons.masses[permutation], phonons.masses, atol=1e-8, rtol=0):
                raise ValueError("Symmetry must preserve atomic masses")
            transformed = np.empty_like(vectors)
            transformed[:, permutation] = vectors @ rotation.T
            transformed = transformed.reshape(len(group), -1)
            gamma = basis.conj() @ transformed.T
            leakage = max(leakage, np.linalg.norm(transformed - gamma.T @ basis) / np.sqrt(len(group)))
            matrices.append(gamma)
        if strict and leakage > tolerance:
            raise ValueError(f"Modes {[int(i + 1) for i in group]} are not a closed symmetry subspace "
                             f"(leakage {leakage:.3g} > {tolerance:g}); check phonon/structure alignment and irreps")
        blocks.append(np.asarray(matrices))
        leakages.append(leakage)
    return blocks, np.asarray(leakages)




def mode_covariance_residuals(tensors, operations, groups, representations):
    """Direct RMS-over-operations residual of S R_k S^T = sum_p Gamma_pk R_p."""
    tensors = np.asarray(tensors)
    if tensors.ndim != 3 or tensors.shape[1:] != (3, 3) or not np.isfinite(tensors).all():
        raise ValueError("Finite mode tensors with shape (nmodes, 3, 3) required")
    if not operations:
        raise ValueError("Symmetry operations required")
    if len(groups) != len(representations) or sorted(i for g in groups for i in g) != list(range(len(tensors))):
        raise ValueError("Representations must match a partition of all modes")
    squared = np.zeros(len(tensors))
    lhs_squared, rhs_squared = np.zeros(len(tensors)), np.zeros(len(tensors))
    for group, matrices in zip(groups, representations):
        matrices = np.asarray(matrices)
        if matrices.shape != (len(operations), len(group), len(group)) or not np.isfinite(matrices).all():
            raise ValueError("Invalid partner representation matrices")
        for (_, rotation), gamma in zip(operations, matrices):
            lhs = rotation @ tensors[group] @ rotation.T
            rhs = np.einsum("pk,pij->kij", gamma, tensors[group])
            squared[group] += np.sum(np.abs(lhs - rhs) ** 2, axis=(1, 2)) / len(operations)
            lhs_squared[group] += np.sum(np.abs(lhs) ** 2, axis=(1, 2)) / len(operations)
            rhs_squared[group] += np.sum(np.abs(rhs) ** 2, axis=(1, 2)) / len(operations)
    def result(a, b, c):
        return dict(absolute_residual=np.sqrt(a), squared_residual=a,
                    comparison_norm=np.maximum(np.sqrt(b), np.sqrt(c)),
                    lhs_norm=np.sqrt(b), rhs_norm=np.sqrt(c))
    per_mode = result(squared, lhs_squared, rhs_squared)
    # Summing partner squared norms gives a basis-independent multiplet diagnostic.
    multiplets = result(np.array([squared[g].sum() for g in groups]),
                        np.array([lhs_squared[g].sum() for g in groups]),
                        np.array([rhs_squared[g].sum() for g in groups]))
    return per_mode, multiplets


def angular_intensities(tensors, angles_degrees, surface_normal="z", scan="parallel"):
    """Return |e_s^dagger R e_i|^2 for rotating linear polarizations (nmodes, nangles)."""
    if surface_normal not in ("x", "y", "z") or scan not in ("parallel", "crossed"):
        raise ValueError("surface_normal must be x/y/z and scan parallel/crossed")
    tensors = np.asarray(tensors)
    if tensors.ndim == 2:
        tensors = tensors[None]
    if tensors.ndim != 3 or tensors.shape[1:] != (3, 3) or not np.isfinite(tensors).all():
        raise ValueError("Tensors must be finite 3x3 matrices")
    angles = np.deg2rad(np.asarray(angles_degrees, dtype=float))
    if angles.ndim != 1 or not np.isfinite(angles).all():
        raise ValueError("Angles must be a finite one-dimensional array")
    first, second = {"x": (1, 2), "y": (2, 0), "z": (0, 1)}[surface_normal]
    incident = np.zeros((len(angles), 3))
    incident[:, first], incident[:, second] = np.cos(angles), np.sin(angles)
    if scan == "parallel":
        scattered = incident
    else:
        scattered = np.zeros_like(incident)
        scattered[:, first], scattered[:, second] = -np.sin(angles), np.cos(angles)
    return np.abs(np.einsum("ti,mij,tj->mt", scattered, tensors, incident)) ** 2


def _weights(frequencies, cutoff):
    weights = np.ones_like(frequencies)
    if cutoff is not None:
        selected = (frequencies > 0) & (frequencies >= cutoff / CM_PER_THZ)
        occupation = 1 / np.expm1(frequencies[selected] * CM_PER_THZ * 0.004824125)
        weights[selected] = (occupation + 1) / frequencies[selected]
    return weights


def _write_csv(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def _plot_style(style):
    import matplotlib.pyplot as plt
    import shutil
    if style is not None:
        plt.style.use(style)
    if shutil.which("latex") is None:
        plt.rcParams["text.usetex"] = False
    return plt


def plot_absolute_residual(ax, frequencies, residuals, intensities, color='C0'):
    """Intensity-sized opaque dots with short horizontal ticks."""
    from matplotlib.colors import to_rgba
    x, y, strength = (np.asarray(v, dtype=float) for v in (frequencies, residuals, intensities))
    if x.shape != y.shape or x.shape != strength.shape or x.ndim != 1:
        raise ValueError('Frequency, residual and intensity arrays must match')
    if not np.isfinite(strength).all() or np.any(strength < 0):
        raise ValueError('Intensity weights must be finite and nonnegative')
    scale = strength.max(initial=0)
    weight = np.sqrt(strength / scale) if scale else np.zeros_like(strength)
    rgba = np.tile(to_rgba(color), (len(x), 1))
    rgba[:, 3] = 1.
    valid = np.isfinite(x) & np.isfinite(y)
    dots = ax.scatter(x[valid], y[valid], s=(5 + 45 * weight)[valid], c=rgba[valid], linewidths=0, zorder=3)
    ax.hlines(y[valid], x[valid] - 3, x[valid] + 3, colors=rgba[valid], linewidths=.7, zorder=2)
    return dots


def _plot_check(output, energy, frequencies, residuals, style, prefix="symmetry_residuals", intensities=None):
    plt = _plot_style(style)
    fig, ax = plt.subplots(figsize=(14, 6), layout="constrained")
    if intensities is None:
        intensities = np.ones(len(frequencies))
    plot_absolute_residual(ax, frequencies, residuals['absolute_residual'], intensities)
    valid_residuals = np.asarray(residuals['absolute_residual'], dtype=float)
    valid_residuals = valid_residuals[np.isfinite(valid_residuals)]
    if valid_residuals.size:
        mean = valid_residuals.mean()
        ax.axhline(mean, color='black', linestyle='--', linewidth=1.2,
                   label=f'Arithmetic mean = {mean:.4g}')
        ax.legend()
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.set_ylim(bottom=0.)
    ax.set_xlabel(r"Phonon frequency (cm$^{-1}$)")
    ax.set_ylabel('Absolute residual (tensor units)')
    missing = np.count_nonzero(~np.isfinite(residuals["absolute_residual"]))
    note = f"; {missing} invalid subspaces omitted" if missing else ""
    fig.suptitle(f"Laser energy: {energy:.2f} eV{note}")
    fig.savefig(output / f"{prefix}_{energy:.2f}eV.png", dpi=250)
    plt.close(fig)


def _plot_angular(output, energy, frequencies, curves, labels, style):
    plt = _plot_style(style)
    for mode, (angles, values) in curves.items():
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), subplot_kw={"projection": "polar"}, layout="constrained")
        maximum = max(values["original_parallel"].max(), values["original_crossed"].max())
        for ax, scan in zip(axes, ("parallel", "crossed")):
            curve = values[f"original_{scan}"]
            ax.plot(np.deg2rad(angles), curve, color="#2470a0")
            ax.set_ylim(0, 1.15 * maximum if maximum else 1)
            ax.set_title(scan.capitalize())
            ax.tick_params(labelsize=9)
        fig.suptitle(f"Mode {mode}, {labels[mode - 1]}, {frequencies[mode - 1]:.2f} cm⁻¹; {energy:.2f} eV\nRadius: Raman intensity (arb. units)")
        fig.savefig(output / f"angular_intensity_mode{mode}_{energy:.2f}eV.png", dpi=250)
        plt.close(fig)


def _mode_analysis(energies, symmetry_average=True):
    settings = read_settings()
    phonons = read_registered_modes(provenance="atom_registration.json")
    nmodes = len(phonons.frequencies_thz)
    if (nmodes != 3 * phonons.natoms or not np.isfinite(phonons.frequencies_thz).all()
            or not np.isfinite(phonons.eigendisplacements).all() or np.any(phonons.masses <= 0)):
        raise ValueError("Invalid Gamma frequencies, eigenvectors, or atomic masses")
    if symmetry_average:
        structure_path = "ref_poscar.vasp" if Path("ref_poscar.vasp").is_file() else "CONTCAR"
        structure = read_structure(structure_path)
        if structure.natoms != phonons.natoms:
            raise ValueError("Structure and band.yaml atom counts differ")
        operations = symmetry_operations(structure)
    labels = read_irreps() or []
    labels = (labels + ["---"] * nmodes)[:nmodes]
    energies = settings.laser_energies if energies is None else energies
    if not energies or not np.isfinite(energies).all():
        raise ValueError("At least one finite laser energy is required")
    for energy in energies:
        real, imaginary = read_dielectric_derivatives(f"epsilon_derivative_{energy:.2f}", phonons.natoms)
        derivatives = real + 1j * imaginary
        original_modes = mode_tensors(derivatives, phonons.eigendisplacements)
        projected_modes = None
        if symmetry_average:
            projected = project_derivatives(derivatives, operations)
            if np.linalg.norm(project_derivatives(projected, operations) - projected) > 1e-5 * max(np.linalg.norm(projected), 1e-30):
                raise ValueError("Group average is not idempotent; check symmetry operations and structure")
            projected_modes = mode_tensors(projected, phonons.eigendisplacements)
        yield energy, phonons, original_modes, projected_modes, labels, settings


def run_check(energies=None, output="check", epsilon_fraction=1e-12, plot=False, style=None):
    """Export all mode symmetry residuals without altering the supplied derivatives."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for energy, phonons, original_modes, projected_modes, labels, _ in _mode_analysis(energies):
        frequencies = phonons.frequencies_thz * CM_PER_THZ
        residuals = projection_residuals(original_modes, projected_modes, epsilon_fraction)
        keys = ("absolute_residual", "squared_residual", "relative_residual", "original_norm", "projected_norm")
        _write_csv(output / f"symmetry_residuals_{energy:.2f}eV.csv",
                   ("mode", "frequency_cm-1", "irrep", *keys, "epsilon"),
                   ((index + 1, frequency, labels[index], *(residuals[key][index] for key in keys), residuals["epsilon"])
                    for index, frequency in enumerate(frequencies)))
        if plot:
            _plot_check(output, energy, frequencies, residuals, style)
        print(f"Wrote symmetry checks for {energy:.2f} eV to {output}")
    (output / "README.txt").write_text(
        "Symmetry consistency of the supplied final atomic derivatives.\n"
        "Projection averages all structural symmetry operations before mode contraction.\n"
        "Residual = norm(original - projected); relative denominator = max(norms) + epsilon.\n"
        "epsilon = epsilon_fraction * max(1, largest original/projected mode norm) per laser energy.\n"
        "Atomic derivatives already include the legacy volume/(4*pi) factor; it is not applied again.\n"
        "Reconstruction can enforce some constraints; this diagnostic is not total numerical error.\n"
        "Electronic/magnetic symmetry compatibility is assumed, not tested.\n"
        "Mode eigenvectors are held fixed; atomic projection alone does not enforce mode symmetry.\n"
    )


def _prepare_mode_check(subspace_tolerance=1e-5):
    """Prepare frequency-independent data; callers must validate files before reuse."""
    phonons = read_gamma_modes(complex_eigenvectors=True)
    structure = read_structure("ref_poscar.vasp" if Path("ref_poscar.vasp").is_file() else "CONTCAR")
    if structure.natoms != phonons.natoms:
        raise ValueError("Structure and band.yaml atom counts differ")
    phonons, atom_mapping, origin_shift, alignment_error = align_mode_atoms(phonons, structure)
    operations = symmetry_operations(structure)
    with Path("irreps.yaml").open() as stream:
        irreps = yaml.load(stream, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    if not np.allclose(irreps["q-position"], 0, atol=1e-8):
        raise ValueError("irreps.yaml must describe Gamma modes")
    entries = irreps["normal_modes"]
    groups = [np.asarray(entry["band_indices"], dtype=int) - 1 for entry in entries]
    representations, leakages = mode_representations(phonons, operations, groups, subspace_tolerance, strict=False)
    statuses = ["ok" if leakage <= subspace_tolerance else "nonclosed_subspace" for leakage in leakages]
    # Validate characters when supplied, matching rotations rather than relying on order.
    for group_index, entry in enumerate(entries):
        if "characters" not in entry:
            continue
        if len(entry["characters"]) != len(irreps["rotations"]):
            raise ValueError("Irrep character count differs from rotation count")
        for rotation_entry, character in zip(irreps["rotations"], entry["characters"]):
            cartesian = structure.lattice.T @ np.asarray(rotation_entry["matrix"]) @ np.linalg.inv(structure.lattice.T)
            matches = [i for i, (_, rotation) in enumerate(operations) if np.allclose(rotation, cartesian, atol=1e-6)]
            expected = character[0] * np.exp(1j * np.deg2rad(character[1]))
            if not matches or any(abs(np.trace(representations[group_index][i]) - expected) > 1e-3 for i in matches):
                statuses[group_index] = "character_mismatch" if statuses[group_index] == "ok" else statuses[group_index]
    mode_statuses = [None] * len(phonons.frequencies_thz)
    for group, status in zip(groups, statuses):
        for index in group:
            mode_statuses[index] = status
    invalid = [j for j, status in enumerate(statuses) if status != "ok"]
    if invalid:
        print(f"Warning: {len(invalid)} partner groups fail phonon-subspace/character validation; "
              "their tensor residuals are unavailable (NaN), not plotted.")
    if len(invalid) == len(groups):
        raise ValueError("No validated mode subspaces; refusing to create misleading residual plots")
    return (phonons, atom_mapping, origin_shift, alignment_error, operations, entries,
            groups, representations, leakages, statuses, mode_statuses, invalid)


def run_mode_check(energies=None, output="check", plot=False,
                   style=None, subspace_tolerance=1e-5, _prepared=None):
    """Check mode covariance directly, independently of atomic projection.

    The private prepared-data hook is for batch callers that first verify
    identical band.yaml, irreps.yaml, symmetry, and reference structure files.
    """
    (phonons, atom_mapping, origin_shift, alignment_error, operations, entries,
     groups, representations, leakages, statuses, mode_statuses, invalid) = (
        _prepare_mode_check(subspace_tolerance) if _prepared is None else _prepared)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    keys = ("absolute_residual", "squared_residual", "lhs_norm", "rhs_norm", "comparison_norm")
    _write_csv(output / "mode_check_atom_mapping.csv", ("phonon_atom", "derivative_atom"),
               ((i + 1, int(target + 1)) for i, target in enumerate(atom_mapping)))
    settings = read_settings()
    energies = settings.laser_energies if energies is None else energies
    if not energies or not np.isfinite(energies).all():
        raise ValueError("At least one finite laser energy is required")
    labels = ["---"] * len(phonons.frequencies_thz)
    for entry, group in zip(entries, groups):
        for index in group:
            labels[index] = entry.get("ir_label", "---")
    for energy in energies:
        real, imaginary = read_dielectric_derivatives(f"epsilon_derivative_{energy:.2f}", phonons.natoms)
        tensors = mode_tensors(real + 1j * imaginary, phonons.eigendisplacements)
        residuals, multiplets = mode_covariance_residuals(tensors, operations, groups, representations)
        for j in invalid:
            for key in keys:
                residuals[key][groups[j]] = np.nan
                multiplets[key][j] = np.nan
        frequencies = phonons.frequencies_thz * CM_PER_THZ
        _write_csv(output / f"mode_covariance_residuals_{energy:.2f}eV.csv",
                   ("mode", "frequency_cm-1", "irrep", "status", *keys),
                   ((i + 1, f, labels[i], mode_statuses[i], *(residuals[key][i] for key in keys))
                    for i, f in enumerate(frequencies)))
        group_frequencies = np.array([frequencies[g].mean() for g in groups])
        _write_csv(output / f"multiplet_covariance_residuals_{energy:.2f}eV.csv",
                   ("modes", "frequency_cm-1", "irrep", "status", "subspace_leakage", *keys),
                   ((" ".join(str(i + 1) for i in g), group_frequencies[j], entries[j].get("ir_label", "---"),
                     statuses[j], leakages[j], *(multiplets[key][j] for key in keys))
                    for j, g in enumerate(groups)))
        if plot:
            axes = {'z': (0, 1), 'y': (0, 2), 'x': (1, 2)}.get(settings.surface_normal, (0, 1, 2))
            intensity = sum(np.abs(tensors[:, i, j]) ** 2 for i in axes for j in axes) * _weights(phonons.frequencies_thz, settings.frequency_cutoff_cm1)
            _plot_check(output, energy, frequencies, residuals, style, "mode_covariance_residuals", intensity)
            _plot_check(output, energy, group_frequencies, multiplets, style, "multiplet_covariance_residuals", np.array([intensity[g].sum() for g in groups]))
        print(f"Wrote direct mode-covariance checks for {energy:.2f} eV to {output}")
    (output / "MODE_CHECK_README.txt").write_text(
        "Direct constraint: S_g R_k S_g^T = sum_p Gamma_pk(g) R_p. No atomic projection.\n"
        "Gamma_pk = <e_p|T(g)e_k> uses the band.yaml mass-weighted eigenvector basis.\n"
        "Phonon atoms are aligned to derivative atoms; mapping is saved separately.\n"
        "Spectrum and angular contractions also require geometry-based atom registration.\n"
        f"Fractional origin shift (phonon to derivative): {origin_shift.tolist()}; max mismatch {alignment_error:g} Angstrom.\n"
        "irreps.yaml supplies partner groups; closure and supplied characters are validated.\n"
        f"Subspace tolerance = {subspace_tolerance:g}. Invalid groups have status flags and NaN residuals.\n"
        "Absolute residual = sqrt(mean_g ||lhs_g-rhs_g||_F^2); squared residual is its square.\n"
        "Comparison magnitude = max(RMS lhs norm, RMS rhs norm). No relative normalization or regularizer.\n"
        "Multiplet norms sum squared partner norms; these are invariant under unitary partner mixing.\n"
        "Individual partner residuals depend on the eigenvector basis. Zero tensors have zero residual.\n"
        "Structural symmetries must preserve the electronic/magnetic state; this is assumed.\n"
        "This checks covariance, not total numerical accuracy; derivatives are never changed.\n"
    )


def run_angular(modes, energies=None, output="angular", angle_step=1.0,
                surface_normal=None, plot=False, style=None):
    """Export selected mode polarization scans, independently of run_check."""
    if not np.isfinite(angle_step) or not 0 < angle_step <= 360:
        raise ValueError("angle_step must be between 0 and 360 degrees")
    selected_modes = list(dict.fromkeys(modes))
    if not selected_modes:
        raise ValueError("At least one mode number is required")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    angles = np.append(np.arange(0, 360, angle_step), 360.0)
    for energy, phonons, original_modes, _, labels, settings in _mode_analysis(energies, symmetry_average=False):
        if any(index < 1 or index > len(phonons.frequencies_thz) for index in selected_modes):
            raise ValueError(f"Mode numbers must lie between 1 and {len(phonons.frequencies_thz)}")
        normal = surface_normal or settings.surface_normal or "z"
        if normal not in ("x", "y", "z"):
            raise ValueError("Surface normal must be x, y, or z")
        weights = _weights(phonons.frequencies_thz, settings.frequency_cutoff_cm1)
        frequencies = phonons.frequencies_thz * CM_PER_THZ
        curves = {}
        for mode in selected_modes:
            values = {}
            for scan in ("parallel", "crossed"):
                values[f"original_{scan}"] = angular_intensities(original_modes[mode - 1], angles, normal, scan)[0] * weights[mode - 1]
            keys = tuple(values)
            _write_csv(output / f"angular_intensity_mode{mode}_{energy:.2f}eV.csv",
                       ("angle_degrees", *keys), zip(angles, *(values[key] for key in keys)))
            curves[mode] = (angles, values)
        if plot:
            _plot_angular(output, energy, frequencies, curves, labels, style)
        print(f"Wrote {len(curves)} angular scans for {energy:.2f} eV to {output}")
    (output / "README.txt").write_text(
        f"Angular scans use the {normal} normal, rotating both incident polarization and analyzer.\n"
        "For normals x/y/z, angle zero and positive 90 degrees are y/z, z/x, x/y respectively.\n"
        "Intensity = |e_s^dagger R e_i|^2, with the input cutoff/298 K frequency weighting.\n"
        "No angular normalization or irrep filtering. Radii are in the existing intensity units.\n"
        "Degenerate modes are exported individually; sum their intensities for an unresolved multiplet.\n"
        "Scans use the calculated complex mode tensors without symmetry projection.\n"
    )
