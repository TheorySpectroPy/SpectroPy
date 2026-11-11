"""Numerical symmetry operations and partner representations for Raman calculations."""

from pathlib import Path
import numpy as np
import yaml

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
