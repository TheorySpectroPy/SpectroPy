"""Numerical symmetry operations and partner representations for Raman calculations."""

import numpy as np
from .symmetry import symmetry_operations


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
