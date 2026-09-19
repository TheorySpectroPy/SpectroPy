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
    by_axis = {}
    for direction, epsilon in measured:
        axis, sign = _classify_axis(direction / np.linalg.norm(direction), atol)
        by_axis.setdefault(axis, {})[sign] = (epsilon, np.linalg.norm(direction))
    slices = []
    for axis, sides in by_axis.items():
        direction = np.eye(3)[axis]
        if 1 in sides and -1 in sides:
            positive, negative = sides[1][0], sides[-1][0]
        elif 1 in sides:
            positive = sides[1][0]
            negative = synthesize_missing_sign(positive, direction, site_symmetry_cart, atol)
        else:
            negative = sides[-1][0]
            positive = synthesize_missing_sign(negative, -direction, site_symmetry_cart, atol)
        # Every finite difference is divided by its own step length, so the
        # derivatives do not depend on the spread of the microscopic step sizes.
        step = np.mean([magnitude for _, magnitude in sides.values()])
        slices.append((direction, (positive - negative) / (2 * step)))
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


# ---------------------------------------------------------------------------
# Fractional-basis mode ("fractional") support
#
# The `fractional` displacement mode probes along lattice-vector directions
# (a, b, a+c, ...) instead of Cartesian axes.  VASP still returns the
# dielectric tensor in Cartesian coordinates, so each central-difference
# measurement yields
#
#     S_f = (eps(+d) - eps(-d)) / 2  ~  D_cart . d  =  D_frac . f
#
# where f is the fractional displacement direction, d = f @ lattice is the
# Cartesian displacement, and the per-atom derivative tensors are
#
#     D_cart[a,b,g] = d eps_ab / d r_g    (displacement index in Cartesian)
#     D_frac[a,b,c] = d eps_ab / d u_c    (displacement index in fractional)
#
# with r = u @ lattice.  The measured slice S_f is therefore the derivative
# contracted with a *fractional* displacement index.
# ---------------------------------------------------------------------------


def reconstruct_fractional_tensor(measured, site_symmetry_frac, lattice):
    """Reconstruct the per-atom Cartesian derivative tensor from fractional directions.

    PSEUDOCODE (implementation intentionally left blank for now):

    measured : list of (fractional_direction f, slice S_f), where
        S_f = (eps(+d) - eps(-d)) / 2   (Cartesian 3x3 matrix),
        with d = f @ lattice the Cartesian displacement.

    1. Expand the orbit of every measured direction under the site symmetry
       using the MIXED transformation: the displacement index rotates with the
       FRACTIONAL rotation R, while the two dielectric indices rotate with the
       CARTESIAN rotation R_cart = frac_to_cart_rotation(R, lattice):

           for R in site_symmetry_frac:
               f'   = f @ R.T                       # fractional direction image
               S_f' = R_cart @ S_f @ R_cart.T       # eps rotates in Cartesian

       (This is the key difference from `reconstruct_full_tensor`, which uses
       the same Cartesian R for both the direction and the eps indices.)

    2. Collect the orbit pairs (f', S_f') into a dict keyed by direction.

    3. Select 3 linearly independent fractional directions {f1, f2, f3} from
       the orbit (basis selection as in `reconstruct_full_tensor`).

    4. Assemble D_frac (3x3x3) by expressing each fractional axis as a linear
       combination of the basis directions, exactly as `reconstruct_full_tensor`
       assembles D_cart, using the relation D_frac . f_k = S_{f_k}.

    5. Convert to the Cartesian displacement basis:
           D_cart = convert_fractional_tensor_to_cartesian(D_frac, lattice)

    6. Return D_cart (3x3x3), the same shape `compute_atom_tensor` returns.
    """
    raise NotImplementedError("fractional-mode tensor reconstruction is not implemented yet")


def convert_fractional_tensor_to_cartesian(D_frac, lattice):
    """Convert a per-atom derivative tensor from the fractional displacement
    basis to the Cartesian displacement basis.

    PSEUDOCODE (implementation intentionally left blank for now):

    Definitions
    -----------
        D_frac[a, b, c] = d eps_ab / d u_c      (u = fractional displacement)
        D_cart[a, b, g] = d eps_ab / d r_g      (r = Cartesian displacement)
        r = u @ lattice,  i.e.  r_g = sum_c u_c lattice[c, g]

    Chain rule
    ----------
        d eps_ab / d u_c = sum_g (d eps_ab / d r_g) (d r_g / d u_c)
                         = sum_g D_cart[a, b, g] lattice[c, g]

    i.e. for every dielectric index pair (a, b):
        D_frac[a, b, :] = lattice @ D_cart[a, b, :]

    Inverting the 3x3 relation for each (a, b):
        D_cart[a, b, :] = inv(lattice) @ D_frac[a, b, :]

    Implementation sketch
    ---------------------
        inv_lattice = np.linalg.inv(lattice)
        D_cart = np.zeros_like(D_frac)
        for a in range(3):
            for b in range(3):
                D_cart[a, b, :] = inv_lattice @ D_frac[a, b, :]
        return D_cart

    Notes
    -----
    * Only the displacement index (the last index) is transformed; the two
      dielectric indices (a, b) are already Cartesian (VASP convention).
    * ``lattice`` is the row-vector lattice matrix (rows = a, b, c), the same
      convention as ``Structure.lattice`` and ``spectropy_displacements``.
    * Once every atom's D_cart is known, the Raman tensor for a mode is built
      from D_cart and the Cartesian eigenvectors as usual, so no further
      conversion is required at the spectrum stage.
    """
    raise NotImplementedError("fractional -> Cartesian conversion is not implemented yet")
