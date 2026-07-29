"""
Reconstructs the full per-atom dielectric-derivative tensor D_ijk (i,j the
dielectric-tensor indices, k the Cartesian displacement direction) from the
reduced set of directions phonopy's generate_displacements() decided were
sufficient (see generate_minimal_displacements.py), plus expands from
symmetry-independent atoms to all atoms.

This is the piece calculate_spectrum.py's read_dielectric_derivatives()
explicitly flagged as missing:
    "Symmetry was used. Full tensor reconstruction is not yet implemented."

Two distinct symmetry expansions are involved, using two different group
actions -- do not conflate them:

1. Direction reconstruction (this module's core): D_ijk(u) is *linear* in
   the displacement direction u. If an atom's own site symmetry contains an
   operation g mapping a measured direction u_m to some new direction u' =
   g @ u_m, then D(u') = g @ D(u_m) @ g^T is known for free -- no new DFT
   calculation needed. Repeating this over the whole site-symmetry group
   orbit of the measured directions usually spans all of R^3 well before
   using all 3 Cartesian axes (that's exactly why phonopy's
   get_least_displacements can skip directions entirely, e.g. MoS2's Mo
   atom needs only +x and +z -- the y-direction response is recovered here
   via the C3/mirror operations of its site symmetry, not measured).

2. Atom expansion (reuses process_symmetry.py's already-computed
   symmetry_operation_matrices): once an independent atom's full D_ijk is
   known, a *different* symmetry operation (one relating that atom to a
   symmetry-equivalent one, e.g. sigma_h relating MoS2's two S atoms) gives
   the equivalent atom's full D_ijk the same way: D_ijk(atom_b) = R D_ijk(atom_a) R^T R
   (all three indices rotated).
"""
import numpy as np


def rotate_tensor3(R, D):
    """D'_ijk = R_ii' R_jj' R_kk' D_i'j'k' for a rank-3 tensor D (3,3,3)."""
    return np.einsum('ip,jq,kr,pqr->ijk', R, R, R, D)


def rotate_slice(R, S):
    """S'_ij = R_ii' R_jj' S_i'j' for a rank-2 tensor slice S (3,3)."""
    return R @ S @ R.T


def _direction_index(d, known_dirs, atol=1e-6):
    for i, d2 in enumerate(known_dirs):
        if np.allclose(np.cross(d, d2), 0, atol=atol):
            return i
    return None


def expand_direction_orbit(measured, site_symmetry_cart):
    """
    measured: list of (direction (3,) unit Cartesian vector, slice (3,3))
        pairs -- one entry per Cartesian direction that was actually
        computed (e.g. finite-differenced from +u/-u VASP runs, or a
        one-sided derivative already resolved against its own symmetry
        partner -- see resolve_one_sided_directions below).
    site_symmetry_cart: (n, 3, 3) array of the atom's own site-symmetry
        rotations in Cartesian coordinates.

    Returns the same list extended with every (direction, slice) pair
    reachable by applying the site-symmetry group to what's already known
    (a direction parallel or antiparallel to one already present is not
    re-added -- D(-u) = -D(u) is handled directly at basis-solve time).
    """
    known = list(measured)
    known_dirs = [d for d, _ in known]
    frontier = list(known)
    while frontier:
        new_frontier = []
        for R in site_symmetry_cart:
            for d, S in frontier:
                d2 = R @ d
                if _direction_index(d2, known_dirs) is None:
                    S2 = rotate_slice(R, S)
                    known.append((d2, S2))
                    known_dirs.append(d2)
                    new_frontier.append((d2, S2))
        frontier = new_frontier
    return known


def reconstruct_full_tensor(measured, site_symmetry_cart):
    """
    Reconstructs D_ijk (shape (3,3,3), last axis = displacement direction k
    in the lab x,y,z basis) for one atom from a partial set of measured
    directions, using that atom's site symmetry to fill in the rest.

    Raises ValueError if the measured directions plus their site-symmetry
    orbit do not span R^3 -- i.e. the input displacement set was
    insufficient (a bug in the caller, not a numerical accuracy issue).
    """
    known = expand_direction_orbit(measured, site_symmetry_cart)
    dirs = np.array([d for d, _ in known])

    basis_idx = []
    basis_rows = np.zeros((0, 3))
    for i, d in enumerate(dirs):
        trial = np.vstack([basis_rows, d])
        if np.linalg.matrix_rank(trial, tol=1e-6) > basis_rows.shape[0]:
            basis_rows = trial
            basis_idx.append(i)
        if len(basis_idx) == 3:
            break
    if len(basis_idx) < 3:
        raise ValueError(
            "Measured directions and their site-symmetry orbit only span "
            f"a {len(basis_idx)}-dimensional subspace -- cannot reconstruct "
            "the full tensor. This means an insufficient displacement set "
            "was used (a generation bug), not a numerical precision issue."
        )

    B = dirs[basis_idx].T  # columns = the 3 basis directions
    B_inv = np.linalg.inv(B)
    S_basis = [known[i][1] for i in basis_idx]

    D_full = np.zeros((3, 3, 3))
    for k, e_k in enumerate(np.eye(3)):
        c = B_inv @ e_k
        D_full[:, :, k] = sum(c[i] * S_basis[i] for i in range(3))
    return D_full


def expand_to_equivalent_atoms(D_by_independent_atom, atom_mapping, mapping_matrices_cart):
    """
    D_by_independent_atom: {atom_idx (1-based, independent atoms only): D (3,3,3)}
    atom_mapping: {atom_idx (1-based): representative_idx (1-based)}, as
        read by process_symmetry.read_symmetry_file (phonopy's atom_mapping).
    mapping_matrices_cart: {(rep_idx, atom_idx): R (3,3) Cartesian rotation}
        -- the Cartesian versions of what process_symmetry.run_mapping()
        already computes and writes to symmetry_operation_matrices (in
        fractional coordinates there; convert with frac_to_cart_rotation
        before calling this).

    Returns {atom_idx (1-based, ALL atoms): D (3,3,3)}.
    """
    D_all = {}
    for atom_idx, rep_idx in atom_mapping.items():
        D_rep = D_by_independent_atom[rep_idx]
        if atom_idx == rep_idx:
            D_all[atom_idx] = D_rep
        else:
            R = mapping_matrices_cart[(rep_idx, atom_idx)]
            D_all[atom_idx] = rotate_tensor3(R, D_rep)
    return D_all


def frac_to_cart_rotation(R_frac, lattice):
    """Cartesian rotation for fractional-coordinate convention x' = R_frac @ x,
    with `lattice` rows = a, b, c (cart = frac @ lattice) -- same convention
    used throughout this pipeline (see visualize_site_symmetry.py)."""
    L = lattice
    return L.T @ R_frac @ np.linalg.inv(L.T)


def _classify_axis(direction, atol=1e-4):
    """direction: unit Cartesian vector, expected axis-aligned (phonopy's
    generate_displacements(is_diagonal=False) only produces +/-x, +/-y, +/-z).
    Returns (axis 0/1/2, sign +1/-1)."""
    axis = int(np.argmax(np.abs(direction)))
    sign = 1 if direction[axis] > 0 else -1
    if not np.allclose(np.abs(direction), np.eye(3)[axis], atol=atol):
        raise ValueError(f"direction {direction} is not axis-aligned")
    return axis, sign


def synthesize_missing_sign(eps_plus, u_axis, site_symmetry_cart, atol=1e-4):
    """
    Given the ABSOLUTE dielectric tensor eps_plus measured at a displacement
    of +u_axis (a Cartesian basis vector, e.g. e_x), finds a site-symmetry
    operation g with g @ u_axis = -u_axis and returns g @ eps_plus @ g^T --
    the absolute dielectric tensor at -u_axis, with NO new DFT calculation.

    This is valid for the full (not just linear-derivative) epsilon, since
    the equilibrium epsilon(0) term is itself separately invariant under any
    site-symmetry operation, so eps(g u) = g @ eps(u) @ g^T holds for the
    whole tensor including its constant part -- not only for D_ijk itself.
    This is exactly the relation phonopy's is_minus_displacement checked for
    before deciding a minus displacement was unnecessary.

    Raises ValueError if no such operation exists in the given site symmetry
    (would mean phonopy's own decision to skip -u_axis was inconsistent with
    this atom's site symmetry -- a bug upstream, not expected in practice).
    """
    for R in site_symmetry_cart:
        if np.allclose(R @ u_axis, -u_axis, atol=atol):
            return rotate_slice(R, eps_plus)
    raise ValueError(
        f"No site-symmetry operation maps {u_axis} to its negative -- "
        "cannot synthesize the missing displacement side. This means the "
        "generator's decision to skip it was wrong for this site symmetry."
    )


def compute_atom_tensor(measured_raw, site_symmetry_cart, atol=1e-4):
    """
    measured_raw: list of (direction (3,) unit Cartesian vector, eps_absolute
        (3,3)) pairs, one per VASP calculation actually run for this atom --
        i.e. raw absolute dielectric tensors, NOT yet finite-differenced.
        Directions must be axis-aligned (see generate_minimal_displacements.py).
    site_symmetry_cart: this atom's site-symmetry rotations, Cartesian.

    Returns D_full (3,3,3): the atom's complete dielectric-derivative tensor,
    with any axis phonopy's generator dropped a +/- pair for reconstructed
    via synthesize_missing_sign, and any axis dropped ENTIRELY reconstructed
    via reconstruct_full_tensor's site-symmetry orbit expansion.
    """
    magnitude = np.mean([np.linalg.norm(d) for d, _ in measured_raw])
    by_axis = {}  # axis -> {sign: eps_absolute}
    for direction, eps in measured_raw:
        unit = direction / np.linalg.norm(direction)
        axis, sign = _classify_axis(unit, atol=atol)
        by_axis.setdefault(axis, {})[sign] = eps

    measured_slices = []
    for axis, sides in by_axis.items():
        e_axis = np.eye(3)[axis]
        if 1 in sides and -1 in sides:
            eps_plus, eps_minus = sides[1], sides[-1]
        elif 1 in sides:
            eps_plus = sides[1]
            eps_minus = synthesize_missing_sign(eps_plus, e_axis, site_symmetry_cart, atol=atol)
        else:
            eps_minus = sides[-1]
            eps_plus = synthesize_missing_sign(eps_minus, -e_axis, site_symmetry_cart, atol=atol)
        slice_axis = (eps_plus - eps_minus) / (2 * magnitude)
        measured_slices.append((e_axis, slice_axis))

    if len(measured_slices) == 3:
        D_full = np.zeros((3, 3, 3))
        for e_axis, S in measured_slices:
            D_full[:, :, int(np.argmax(e_axis))] = S
        return D_full

    return reconstruct_full_tensor(measured_slices, site_symmetry_cart)
