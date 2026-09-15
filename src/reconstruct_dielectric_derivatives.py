import numpy as np


def rotate_tensor3(R, D):
    return np.einsum("ip,jq,kr,pqr->ijk", R, R, R, D)


def rotate_slice(R, S):
    return R @ S @ R.T


def frac_to_cart_rotation(R_frac, lattice):
    return lattice.T @ R_frac @ np.linalg.inv(lattice.T)


def _direction_index(direction, known_directions, atol=1e-6):
    for index, known in enumerate(known_directions):
        if np.allclose(np.cross(direction, known), 0, atol=atol):
            return index
    return None


def expand_direction_orbit(measured, site_symmetry_cart):
    known = list(measured)
    known_directions = [direction for direction, _ in known]
    frontier = list(known)
    while frontier:
        next_frontier = []
        for rotation in site_symmetry_cart:
            for direction, tensor in frontier:
                mapped_direction = rotation @ direction
                if _direction_index(mapped_direction, known_directions) is None:
                    mapped_tensor = rotate_slice(rotation, tensor)
                    known.append((mapped_direction, mapped_tensor))
                    known_directions.append(mapped_direction)
                    next_frontier.append((mapped_direction, mapped_tensor))
        frontier = next_frontier
    return known


def reconstruct_full_tensor(measured, site_symmetry_cart):
    known = expand_direction_orbit(measured, site_symmetry_cart)
    directions = np.array([direction for direction, _ in known])
    basis_indices = []
    basis = np.zeros((0, 3))
    for index, direction in enumerate(directions):
        candidate = np.vstack([basis, direction])
        if np.linalg.matrix_rank(candidate, tol=1e-6) > basis.shape[0]:
            basis = candidate
            basis_indices.append(index)
        if len(basis_indices) == 3:
            break
    if len(basis_indices) < 3:
        raise ValueError("Measured directions and their site-symmetry orbit do not span three dimensions")
    inverse_basis = np.linalg.inv(directions[basis_indices].T)
    slices = [known[index][1] for index in basis_indices]
    result = np.zeros((3, 3, 3))
    for axis, direction in enumerate(np.eye(3)):
        coefficients = inverse_basis @ direction
        result[:, :, axis] = sum(coefficient * tensor for coefficient, tensor in zip(coefficients, slices))
    return result


def _classify_axis(direction, atol=1e-4):
    axis = int(np.argmax(np.abs(direction)))
    sign = 1 if direction[axis] > 0 else -1
    if not np.allclose(np.abs(direction), np.eye(3)[axis], atol=atol):
        raise ValueError(f"direction {direction} is not axis-aligned")
    return axis, sign


def synthesize_missing_sign(epsilon, direction, site_symmetry_cart, atol=1e-4):
    for rotation in site_symmetry_cart:
        if np.allclose(rotation @ direction, -direction, atol=atol):
            return rotate_slice(rotation, epsilon)
    raise ValueError(f"No site-symmetry operation maps {direction} to its negative")


def compute_atom_tensor(measured, site_symmetry_cart, atol=1e-4):
    magnitude = np.mean([np.linalg.norm(direction) for direction, _ in measured])
    by_axis = {}
    for direction, epsilon in measured:
        axis, sign = _classify_axis(direction / np.linalg.norm(direction), atol)
        by_axis.setdefault(axis, {})[sign] = epsilon
    slices = []
    for axis, sides in by_axis.items():
        direction = np.eye(3)[axis]
        if 1 in sides and -1 in sides:
            positive, negative = sides[1], sides[-1]
        elif 1 in sides:
            positive = sides[1]
            negative = synthesize_missing_sign(positive, direction, site_symmetry_cart, atol)
        else:
            negative = sides[-1]
            positive = synthesize_missing_sign(negative, -direction, site_symmetry_cart, atol)
        slices.append((direction, (positive - negative) / (2 * magnitude)))
    if len(slices) == 3:
        result = np.zeros((3, 3, 3))
        for direction, tensor in slices:
            result[:, :, int(np.argmax(direction))] = tensor
        return result
    return reconstruct_full_tensor(slices, site_symmetry_cart)


def expand_to_equivalent_atoms(measured_tensors, atom_mapping, mapping_matrices_cart):
    all_tensors = dict(measured_tensors)
    for atom, representative in atom_mapping.items():
        if atom in all_tensors:
            continue
        tensor = measured_tensors[representative]
        all_tensors[atom] = rotate_tensor3(mapping_matrices_cart[(representative, atom)], tensor)
    return all_tensors


def calculate_derivatives(measurements, lattice, site_symmetries=None):
    site_symmetries = {} if site_symmetries is None else site_symmetries
    real = {}
    imaginary = {}
    for atom, values in measurements.items():
        site_symmetry = site_symmetries.get(atom, {"rotations": []})
        rotations = np.array([frac_to_cart_rotation(rotation, lattice) for rotation in site_symmetry["rotations"]])
        real[atom] = compute_atom_tensor([(direction, tensor) for direction, tensor, _ in values], rotations)
        imaginary[atom] = compute_atom_tensor([(direction, tensor) for direction, _, tensor in values], rotations)
    return real, imaginary


def reconstruct_derivatives(measurements, lattice, site_symmetries, atom_mapping, mapping_matrices_fractional):
    mapping_matrices_cart = {pair: frac_to_cart_rotation(rotation, lattice) for pair, rotation in mapping_matrices_fractional.items()}
    real, imaginary = calculate_derivatives(measurements, lattice, site_symmetries)
    return (
        expand_to_equivalent_atoms(real, atom_mapping, mapping_matrices_cart),
        expand_to_equivalent_atoms(imaginary, atom_mapping, mapping_matrices_cart),
    )
