import json
from pathlib import Path

import numpy as np


def rotate_tensor3(R, D):
    return np.einsum("ip,jq,kr,pqr->ijk", R, R, R, D)


def rotate_slice(R, S):
    return R @ S @ R.T


def frac_to_cart_rotation(R_frac, lattice):
    return lattice.T @ R_frac @ np.linalg.inv(lattice.T)


def project_site_tensor(tensor, rotations):
    """Orthogonal site-group projection of all three Cartesian indices.

    Accepts either consistent rank-three index convention used by SpectroPy;
    the same polar-vector rotation acts on all indices. No reciprocity or ASR
    is added. The caller supplies the complete electronic-state site group.
    """
    tensor = np.asarray(tensor)
    rotations = np.asarray(rotations, dtype=float)
    if tensor.shape != (3, 3, 3) or not np.isfinite(tensor).all():
        raise ValueError('Site projection requires a finite 3x3x3 tensor')
    if rotations.ndim != 3 or rotations.shape[1:] != (3, 3) or len(rotations) == 0:
        raise ValueError('Site projection requires a nonempty site group')
    if not np.allclose(rotations @ rotations.transpose(0, 2, 1), np.eye(3), atol=1e-8, rtol=0):
        raise ValueError('Site rotations must be Cartesian orthogonal matrices')
    projected = sum(rotate_tensor3(r, tensor) for r in rotations) / len(rotations)
    if not all(np.allclose(rotate_tensor3(r, projected), projected, atol=1e-8, rtol=1e-8)
               for r in rotations):
        raise ValueError('Site group does not produce an invariant projection')
    return projected


def project_reciprocal_site_tensor(tensor, rotations):
    """Project D[alpha,beta,i] onto site symmetry and dielectric reciprocity.

    Transpose, not conjugate transpose: real and imaginary parts are each
    symmetric for reciprocal nonmagnetic response. The projectors commute.
    Not applicable to magneto-optical response. No acoustic sum rule is added.
    """
    projected = project_site_tensor(tensor, rotations)
    return .5 * (projected + projected.swapaxes(0, 1))


def project_crystal_tensor_field(tensors, lattice, positions, rotations, translations,
                                 species=None, tolerance=1e-4):
    """Equal-weight least-squares projection of a complete Cartesian field.

    D[atom,alpha,beta,i] includes every atom, with all three Cartesian
    derivative slices. Complete space-group averaging enforces equivalent-atom
    covariance and site constraints simultaneously. Dielectric transpose
    symmetrization adds reciprocity; real/imaginary parts are treated linearly.
    No acoustic sum rule, mode projection, or volume factor is introduced.
    For nonorthogonal/mixed-direction observations use a constrained fit in
    observation space instead, rather than this Cartesian Frobenius loss.
    """
    tensors=np.asarray(tensors);lattice=np.asarray(lattice,dtype=float)
    positions=np.asarray(positions,dtype=float);rotations=np.asarray(rotations)
    translations=np.asarray(translations,dtype=float)
    n=len(positions)
    if tensors.shape!=(n,3,3,3) or not np.isfinite(tensors).all():
        raise ValueError('Crystal projection requires finite Cartesian tensors for every atom')
    if rotations.ndim!=3 or rotations.shape[1:]!=(3,3) or not len(rotations) or translations.shape!=(len(rotations),3):
        raise ValueError('Invalid space-group operations')
    result=np.zeros_like(tensors,dtype=np.result_type(tensors,float))
    for w,t in zip(rotations,translations):
        r=frac_to_cart_rotation(w,lattice)
        if not np.allclose(r@r.T,np.eye(3),atol=1e-8,rtol=0):
            raise ValueError('Space-group Cartesian rotation must be orthogonal')
        delta=(positions@w.T+t)[:,None,:]-positions[None,:,:]
        delta-=np.rint(delta)
        distance=np.linalg.norm(delta@lattice,axis=2)
        if species is not None:
            symbols=np.asarray(species)
            if len(symbols)!=n:raise ValueError('Species length differs from tensor field')
            distance=np.where(symbols[:,None]==symbols[None,:],distance,np.inf)
        permutation=distance.argmin(axis=1)
        if len(set(permutation))!=n or np.any(distance[np.arange(n),permutation]>tolerance):
            raise ValueError('Space-group operation does not bijectively map the atomic structure')
        transformed=np.einsum('ac,bd,ij,ncdj->nabi',r,r,r,tensors)
        result[permutation]+=transformed/len(rotations)
    return .5*(result+result.swapaxes(1,2))


def _direction_index(direction, known_directions, atol=1e-6):
    for index, known in enumerate(known_directions):
        if np.allclose(np.cross(direction, known), 0, atol=atol):
            return index
    return None


def expand_direction_orbit(measured, site_symmetry_cart, direction_rotations=None, transform=rotate_slice):
    # The displacement direction normally rotates with the same Cartesian
    # operation as the dielectric indices; fractional mode passes the
    # fractional-basis rotations here instead.
    if direction_rotations is None:
        direction_rotations = site_symmetry_cart
    known = list(measured)
    known_directions = [direction for direction, _ in known]
    frontier = list(known)
    while frontier:
        next_frontier = []
        for rotation, direction_rotation in zip(site_symmetry_cart, direction_rotations):
            for direction, tensor in frontier:
                mapped_direction = direction_rotation @ direction
                if _direction_index(mapped_direction, known_directions) is None:
                    mapped_tensor = transform(rotation, tensor)
                    known.append((mapped_direction, mapped_tensor))
                    known_directions.append(mapped_direction)
                    next_frontier.append((mapped_direction, mapped_tensor))
        frontier = next_frontier
    return known


def _basis_indices(directions):
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
    return basis_indices


def reconstruct_full_tensor(measured, site_symmetry_cart, direction_rotations=None):
    known = expand_direction_orbit(measured, site_symmetry_cart, direction_rotations)
    directions = np.array([direction for direction, _ in known])
    basis_indices = _basis_indices(directions)
    inverse_basis = np.linalg.inv(directions[basis_indices].T)
    slices = [known[index][1] for index in basis_indices]
    result = np.zeros((3, 3, 3))
    for axis, direction in enumerate(np.eye(3)):
        coefficients = inverse_basis @ direction
        result[:, :, axis] = sum(coefficient * tensor for coefficient, tensor in zip(coefficients, slices))
    return result


def _classify_axis(direction, atol=1e-4):
    direction = np.asarray(direction, dtype=float)
    norm = np.linalg.norm(direction)
    if norm == 0:
        raise ValueError("zero displacement direction")
    unit = direction / norm
    axis = int(np.argmax(np.abs(unit)))
    sign = 1 if unit[axis] > 0 else -1
    transverse = np.delete(np.abs(unit), axis)
    if np.max(transverse) > atol:
        raise ValueError(f"direction {direction} is not axis-aligned")
    return axis, sign


def synthesize_missing_sign(epsilon, direction, site_symmetry_cart, atol=1e-4):
    for rotation in site_symmetry_cart:
        if np.allclose(rotation @ direction, -direction, atol=atol):
            return rotate_slice(rotation, epsilon)
    raise ValueError(f"No site-symmetry operation maps {direction} to its negative")


def _validate_difference(scheme, equilibrium):
    if scheme not in {"central", "forward"}:
        raise ValueError("difference scheme must be central or forward")
    if scheme == "forward" and (equilibrium is None or np.shape(equilibrium) != (3, 3)
                                or not np.isfinite(equilibrium).all()):
        raise ValueError("Forward differences require a finite equilibrium dielectric tensor")


def fit_allowed_tensor(slices, rotations):
    """Fit ONLY measured directional responses in the reciprocal site subspace.

    Unmeasured directions never enter as independent observations. The fitted
    representative is constrained before any equivalent-atom transfer.
    """
    eye = np.eye(27).reshape(27,3,3,3)
    projector = np.column_stack([project_reciprocal_site_tensor(e, rotations).ravel() for e in eye])
    values, vectors = np.linalg.eigh((projector+projector.T)/2)
    basis = vectors[:, values > .5].reshape(3,3,3,-1)
    if basis.shape[-1] == 0:
        return np.zeros((3,3,3), dtype=np.result_type(*[s for _,s in slices]))
    design=[]; target=[]
    for direction, response in slices:
        length=np.linalg.norm(direction)
        response=.5*(response+response.T)  # reciprocal measured response first
        design.append(np.einsum('abip,i->abp',basis,direction/length).reshape(9,-1))
        target.append((response/length).ravel())
    design=np.vstack(design);target=np.concatenate(target)
    coefficients,_,rank,_=np.linalg.lstsq(design,target,rcond=1e-10)
    if rank != basis.shape[-1]:
        raise ValueError('Measured seed directions do not determine the site-allowed tensor')
    return np.einsum('abip,p->abi',basis,coefficients)


def compute_atom_tensor(measured, site_symmetry_cart, atol=1e-4,
                        difference_scheme="central", equilibrium=None, constrained=False):
    _validate_difference(difference_scheme, equilibrium)
    by_axis = {}
    for direction, epsilon in measured:
        axis, sign = _classify_axis(direction / np.linalg.norm(direction), atol)
        by_axis.setdefault(axis, {})[sign] = (epsilon, np.linalg.norm(direction))
    slices = []
    for axis, sides in by_axis.items():
        direction = np.eye(3)[axis]
        if difference_scheme == "forward":
            sign = 1 if 1 in sides else -1
            epsilon, step = sides[sign]
            slices.append((direction, (epsilon - equilibrium) / (sign * step)))
            continue
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
    if constrained:
        return fit_allowed_tensor(slices, site_symmetry_cart)
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


def calculate_derivatives(measurements, lattice, site_symmetries=None,
                          difference_scheme="central", equilibrium=None, constrained=False):
    site_symmetries = {} if site_symmetries is None else site_symmetries
    real = {}
    imaginary = {}
    for atom, values in measurements.items():
        site_symmetry = site_symmetries.get(atom, {"rotations": []})
        rotations = np.array([frac_to_cart_rotation(rotation, lattice) for rotation in site_symmetry["rotations"]])
        real[atom] = compute_atom_tensor([(direction, tensor) for direction, tensor, _ in values], rotations,
            difference_scheme=difference_scheme, equilibrium=None if equilibrium is None else equilibrium[0], constrained=constrained)
        imaginary[atom] = compute_atom_tensor([(direction, tensor) for direction, _, tensor in values], rotations,
            difference_scheme=difference_scheme, equilibrium=None if equilibrium is None else equilibrium[1], constrained=constrained)
    return real, imaginary


def reconstruct_derivatives(measurements, lattice, site_symmetries, atom_mapping, mapping_matrices_fractional,
                            difference_scheme="central", equilibrium=None, site_symmetry_projection=False):
    mapping_matrices_cart = {pair: frac_to_cart_rotation(rotation, lattice) for pair, rotation in mapping_matrices_fractional.items()}
    real, imaginary = calculate_derivatives(measurements, lattice, site_symmetries,
                                          difference_scheme, equilibrium, constrained=site_symmetry_projection)
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


def reconstruct_fractional_tensor(measured, site_symmetry_frac, lattice, basis=None):
    """Reconstruct the per-atom Cartesian derivative tensor from fractional directions.

    ``measured`` is a list of (unit fractional direction f, slice D_frac . f).
    A saved basis supplies three response recipes without an orbit search.
    Without one, the legacy mixed-transformation orbit expansion is used.
    D_frac is assembled from three independent directions and its displacement
    index is then converted to Cartesian.
    """
    site_symmetry_cart = [frac_to_cart_rotation(rotation, lattice) for rotation in site_symmetry_frac]
    if basis is None:
        D_frac = reconstruct_full_tensor(measured, site_symmetry_cart, direction_rotations=list(site_symmetry_frac))
    else:
        if not np.allclose([direction for direction, _ in measured], basis["seed_directions"], rtol=0, atol=1e-12):
            raise ValueError("Saved fractional basis does not match the measured directions")
        slices = []
        for entry in basis["entries"]:
            tensor = measured[entry["seed"]][1]
            # Preserve the legacy sequence of rotations, including roundoff.
            for rotation in entry["rotations"]:
                tensor = rotate_slice(np.asarray(rotation), tensor)
            slices.append(tensor)
        inverse_basis = np.asarray(basis["inverse"])
        D_frac = np.zeros((3, 3, 3))
        for axis, direction in enumerate(np.eye(3)):
            coefficients = inverse_basis @ direction
            D_frac[:, :, axis] = sum(coefficient * tensor for coefficient, tensor in zip(coefficients, slices))
    return convert_fractional_tensor_to_cartesian(D_frac, lattice)


def convert_fractional_tensor_to_cartesian(D_frac, lattice):
    """Convert the displacement index from fractional to Cartesian.

    With r = u @ lattice, d eps_ab / d u_c = sum_g D_cart[a, b, g] lattice[c, g],
    so D_cart[a, b, :] = inv(lattice) @ D_frac[a, b, :]. The dielectric indices
    are already Cartesian (VASP convention), as are Phonopy's eigenvectors.
    """
    return np.einsum("gc,abc->abg", np.linalg.inv(lattice), D_frac)


def _fractional_pairs(measured, lattice, atol=1e-6):
    """Group signed measurements in their original order."""
    inverse_lattice = np.linalg.inv(lattice)
    pairs = []
    for direction, epsilon in measured:
        vector = direction @ inverse_lattice
        length = np.linalg.norm(vector)
        if length == 0:
            raise ValueError("Zero fractional displacement")
        unit = vector / length
        for pair in pairs:
            if np.allclose(unit, pair[0], atol=atol):
                pair[2][1] = epsilon
                break
            if np.allclose(unit, -pair[0], atol=atol):
                pair[2][-1] = epsilon
                break
        else:
            pairs.append([unit, length, {1: epsilon}])
    return pairs


def _basis_context(displacements, lattice, site_symmetries):
    atoms = sorted({displacement.atom_index for displacement in displacements})
    return {
        "lattice": np.asarray(lattice).tolist(),
        "displacements": [[d.atom_index, *d.fractional_vector.tolist()] for d in displacements],
        "site_rotations": {str(atom): np.asarray(site_symmetries.get(atom, {"rotations": np.eye(3)[None]})["rotations"]).tolist() for atom in atoms},
    }


def build_fractional_bases(displacements, lattice, site_symmetries):
    """Save the legacy first-independent basis and response recipes, once.

    Use the serialized displacement vectors so generation and reconstruction
    use identical directions despite displacement-file rounding.
    """
    context = _basis_context(displacements, lattice, site_symmetries)
    bases = {}
    for atom, rotations in context["site_rotations"].items():
        measured = [(d.fractional_vector @ lattice, None) for d in displacements if d.atom_index == int(atom)]
        seeds = [pair[0] for pair in _fractional_pairs(measured, lattice)]
        cart = [frac_to_cart_rotation(np.asarray(rotation), lattice) for rotation in rotations]
        recipes = [(seed, {"seed": index, "rotations": []}) for index, seed in enumerate(seeds)]
        def transform(rotation, recipe):
            return {"seed": recipe["seed"], "rotations": recipe["rotations"] + [rotation.tolist()]}
        known = expand_direction_orbit(recipes, cart, [np.asarray(r) for r in rotations], transform)
        directions = np.asarray([direction for direction, _ in known])
        indices = _basis_indices(directions)
        matrix = directions[indices].T
        bases[atom] = {"seed_directions": np.asarray(seeds).tolist(), "matrix": matrix.tolist(),
                       "inverse": np.linalg.inv(matrix).tolist(), "entries": [known[i][1] for i in indices]}
    return {"version": 1, "context": context, "atoms": bases}


def write_fractional_bases(path, bases):
    Path(path).write_text(json.dumps(bases, indent=2) + "\n")


def read_fractional_bases(path, displacements, lattice, site_symmetries):
    bases = json.loads(Path(path).read_text())
    if bases.get("version") != 1 or bases.get("context") != _basis_context(displacements, lattice, site_symmetries):
        raise ValueError("Saved fractional basis is stale; regenerate it for the current displacements, lattice and symmetry")
    for atom in bases["context"]["site_rotations"]:
        basis = bases["atoms"][atom]
        matrix, inverse = np.asarray(basis["matrix"]), np.asarray(basis["inverse"])
        if matrix.shape != (3, 3) or inverse.shape != (3, 3) or not np.allclose(matrix @ inverse, np.eye(3)):
            raise ValueError(f"Invalid saved fractional basis for atom {atom}")
        if len(basis["entries"]) != 3:
            raise ValueError(f"Invalid saved fractional response recipes for atom {atom}")
    return bases


def fractional_finite_differences(measured, site_symmetry_frac, lattice, atol=1e-6,
                                difference_scheme="central", equilibrium=None):
    """Measured seed responses only, before a basis solve or denoising.

    Responses are per unit fractional path parameter. The physical derivative
    is response / norm(unit @ lattice). A signed pair supplies one observation.
    """
    _validate_difference(difference_scheme, equilibrium)
    pairs = _fractional_pairs(measured, lattice, atol)
    slices = []
    for unit, length, sides in pairs:
        generated = False
        if difference_scheme == "forward":
            slices.append((unit, (sides[1] - equilibrium) / length, length, generated))
            continue
        if -1 not in sides:
            for rotation in site_symmetry_frac:
                if np.allclose(rotation @ unit, -unit, atol=atol):
                    sides[-1] = rotate_slice(frac_to_cart_rotation(rotation, lattice), sides[1])
                    generated = True
                    break
            else:
                raise ValueError(f"No site-symmetry operation maps fractional direction {unit} to its negative")
        slices.append((unit, (sides[1] - sides[-1]) / (2 * length), length, generated))
    return slices


def compute_fractional_atom_tensor(measured, site_symmetry_frac, lattice, atol=1e-6, basis=None,
                                   difference_scheme="central", equilibrium=None, constrained=False):
    """Solve measured directional differences in Cartesian tensor components."""
    slices = [(unit, response) for unit, response, _, _ in fractional_finite_differences(
        measured, site_symmetry_frac, lattice, atol, difference_scheme, equilibrium)]
    if constrained:
        return fit_allowed_tensor([(unit @ lattice, response) for unit,response in slices],
            [frac_to_cart_rotation(r,lattice) for r in site_symmetry_frac])
    return reconstruct_fractional_tensor(slices, site_symmetry_frac, lattice, basis=basis)


def reconstruct_fractional_derivatives(measurements, lattice, site_symmetries, atom_mapping, mapping_matrices_fractional, bases=None,
                                       difference_scheme="central", equilibrium=None, site_symmetry_projection=False):
    mapping_matrices_cart = {pair: frac_to_cart_rotation(rotation, lattice) for pair, rotation in mapping_matrices_fractional.items()}
    real = {}
    imaginary = {}
    for atom, values in measurements.items():
        rotations = site_symmetries.get(atom, {"rotations": np.eye(3)[None]})["rotations"]
        basis = None if bases is None else bases["atoms"][str(atom)]
        real[atom] = compute_fractional_atom_tensor([(direction, tensor) for direction, tensor, _ in values], rotations, lattice, basis=basis,
            difference_scheme=difference_scheme, equilibrium=None if equilibrium is None else equilibrium[0], constrained=site_symmetry_projection)
        imaginary[atom] = compute_fractional_atom_tensor([(direction, tensor) for direction, _, tensor in values], rotations, lattice, basis=basis,
            difference_scheme=difference_scheme, equilibrium=None if equilibrium is None else equilibrium[1], constrained=site_symmetry_projection)
    return (
        expand_to_equivalent_atoms(real, atom_mapping, mapping_matrices_cart),
        expand_to_equivalent_atoms(imaginary, atom_mapping, mapping_matrices_cart),
    )
