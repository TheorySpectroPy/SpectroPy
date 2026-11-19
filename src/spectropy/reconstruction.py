import json
from pathlib import Path
from dataclasses import dataclass
from types import SimpleNamespace
from .symmetry import frac_to_cart_rotation, prepare_operations

import numpy as np


def rotate_tensor3(R, D):
    return np.einsum("ip,jq,kr,pqr->ijk", R, R, R, D)


def rotate_slice(R, S):
    return R @ S @ R.T


def project_site_tensor(tensor, rotations):
    """Orthogonal site-group projection of all three Cartesian indices.

    Uses D[alpha,beta,i]; the same polar-vector rotation acts on all indices. No reciprocity or ASR
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
    operations = prepare_operations(lattice, positions, rotations, translations, species, tolerance)
    return project_crystal_field(tensors, operations)


def project_crystal_field(tensors, operations):
    """Project D[atom,alpha,beta,i] using already prepared symmetry operations."""
    tensors = np.asarray(tensors)
    if (tensors.ndim != 4 or tensors.shape[1:] != (3, 3, 3)
            or not np.isfinite(tensors).all()):
        raise ValueError('Crystal projection requires finite Cartesian tensors for every atom')
    result = np.zeros_like(tensors, dtype=np.result_type(tensors, float))
    for permutation, rotation in operations:
        transformed = np.einsum('ac,bd,ij,ncdj->nabi', rotation, rotation, rotation, tensors)
        result[permutation] += transformed / len(operations)
    return .5 * (result + result.swapaxes(1, 2))



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


def _validate_difference(scheme, equilibrium):
    if scheme not in {"central", "forward", "auto"}:
        raise ValueError("difference scheme must be central, forward or auto")
    if scheme == "forward" and (equilibrium is None or np.shape(equilibrium) != (3, 3)
                                or not np.isfinite(equilibrium).all()):
        raise ValueError("Forward differences require a finite equilibrium dielectric tensor")


def allowed_tensor_basis(rotations):
    """Cartesian reciprocal site-allowed basis, prepared once per atom."""
    eye = np.eye(27).reshape(27, 3, 3, 3)
    projector = np.column_stack([project_reciprocal_site_tensor(e, rotations).ravel() for e in eye])
    values, vectors = np.linalg.eigh((projector + projector.T) / 2)
    return vectors[:, values > .5].reshape(3, 3, 3, -1)


def prepare_allowed_fit(directions, rotations):
    """Prepare the existing equal-weight measured-direction least-squares fit."""
    basis = allowed_tensor_basis(rotations)
    lengths = np.linalg.norm(directions, axis=1)
    if basis.shape[-1] == 0:
        return basis, np.zeros((0, len(directions) * 9)), lengths
    design = np.vstack([np.einsum('abip,i->abp', basis, u / length).reshape(9, -1)
                        for u, length in zip(directions, lengths)])
    inverse, _, rank, _ = np.linalg.lstsq(design, np.eye(len(design)), rcond=1e-10)
    if rank != basis.shape[-1]:
        raise ValueError('Measured seed directions do not determine the site-allowed tensor')
    return basis, inverse, lengths


def apply_allowed_fit(responses, fit):
    basis, inverse, lengths = fit
    responses = np.asarray(responses)
    reciprocal = .5 * (responses + responses.swapaxes(1, 2))
    target = (reciprocal / lengths[:, None, None]).ravel()
    return np.einsum('abip,p->abi', basis, inverse @ target)


def fit_allowed_tensor(slices, rotations):
    """Fit ONLY measured directional responses in the reciprocal site subspace.

    Unmeasured directions never enter as independent observations. The fitted
    representative is constrained before any equivalent-atom transfer.
    """
    directions, responses = zip(*slices)
    return apply_allowed_fit(responses, prepare_allowed_fit(directions, rotations))


def _fractional_pairs(measured, lattice, atol=1e-6):
    """Group signed measurements in their original order."""
    inverse_lattice = np.linalg.inv(lattice)
    pairs = []
    for direction, epsilon in measured:
        vector = direction @ inverse_lattice
        length = np.linalg.norm(vector)
        if not np.isfinite(vector).all() or not np.isfinite(length) or length == 0:
            raise ValueError("Directional measurements require finite nonzero displacements; separate equilibrium first")
        unit = vector / length
        for pair in pairs:
            if np.allclose(unit, pair[0], atol=atol, rtol=0):
                raise ValueError("Repeated same-sign displacement direction; use one amplitude at a time")
            if np.allclose(unit, -pair[0], atol=atol, rtol=0):
                if -1 in pair[2]:
                    raise ValueError("Repeated opposite displacement direction")
                if not np.isclose(length, pair[1], atol=1e-10, rtol=1e-8):
                    raise ValueError("Opposite displacements must have equal lengths")
                pair[2][-1] = epsilon
                break
        else:
            pairs.append([unit, length, {1: epsilon}])
    return pairs


def prepare_pairs(directions, site_rotations, lattice, difference_scheme='auto', atol=1e-6,
                  *, equilibrium_available=False):
    """Pair measurement indices and precompute any symmetry-supplied minus sign."""
    if difference_scheme not in {'central', 'forward', 'auto'}:
        raise ValueError('difference scheme must be central, forward or auto')
    pairs = []
    for unit, length, sides in _fractional_pairs(list(zip(directions, range(len(directions)))), lattice, atol):
        negative_rotation = None
        # Equilibrium is a fallback reference, not a reason to discard a
        # symmetry-supplied central pair. Measured minus signs take precedence.
        if difference_scheme in {'central', 'auto'} and -1 not in sides:
            for rotation in site_rotations:
                if np.allclose(rotation @ unit, -unit, atol=atol):
                    negative_rotation = frac_to_cart_rotation(rotation, lattice)
                    break
            else:
                if difference_scheme == 'central':
                    raise ValueError(f'No site-symmetry operation maps fractional direction {unit} to its negative')
        pairs.append((unit, length, sides[1], sides.get(-1), negative_rotation))
    return pairs


def apply_finite_differences(responses, pairs, difference_scheme, equilibrium=None):
    """Complex directional derivatives per unit fractional path parameter."""
    _validate_difference(difference_scheme, equilibrium)
    slices = []
    for unit, length, positive, negative, rotation in pairs:
        forward = difference_scheme == 'forward' or (
            difference_scheme == 'auto' and negative is None and rotation is None)
        if forward:
            _validate_difference('forward', equilibrium)
            derivative = (responses[positive] - equilibrium) / length
        else:
            minus = responses[negative] if negative is not None else rotate_slice(rotation, responses[positive])
            derivative = (responses[positive] - minus) / (2 * length)
        slices.append((unit, derivative, length, rotation is not None))
    return slices


@dataclass(frozen=True)
class AtomicReconstruction:
    """Geometry-only plan; output is complex D[atom,alpha,beta,i] in 1/Angstrom.

    Full-data crystal projection and measured-only constrained fitting remain
    distinct objectives. This plan never changes which observations are fitted.
    """
    natoms: int
    inverse_lattice: np.ndarray
    pairs: dict
    bases: dict
    fits: dict
    atom_mapping: dict
    transfers: dict
    difference_scheme: str

    @property
    def needs_equilibrium(self):
        return self.difference_scheme == 'forward' or (
            self.difference_scheme == 'auto' and any(
                negative is None and rotation is None
                for pairs in self.pairs.values() for _, _, _, negative, rotation in pairs))


def prepare_reconstruction(displacements, lattice, sites, atom_mapping, transfers,
                           difference_scheme='central', denoise=False, bases=None,
                           *, equilibrium_available=False):
    """Build sign mappings, reconstruction recipes, and allowed fits once."""
    if not displacements:
        raise ValueError('No nonzero displacement records to reconstruct')
    pairs, fits = {}, {}
    if not denoise and bases is None:
        bases = build_fractional_bases(displacements, lattice, sites)
    for atom in sorted({d.atom_index for d in displacements}):
        directions = [d.fractional_vector @ lattice for d in displacements if d.atom_index == atom]
        rotations = sites.get(atom, {'rotations': np.eye(3)[None]})['rotations']
        pairs[atom] = prepare_pairs(directions, rotations, lattice, difference_scheme,
                                    equilibrium_available=equilibrium_available)
        if denoise:
            cartesian = [frac_to_cart_rotation(r, lattice) for r in rotations]
            fits[atom] = prepare_allowed_fit([pair[0] @ lattice for pair in pairs[atom]], cartesian)
    return AtomicReconstruction(max(atom_mapping), np.linalg.inv(lattice), pairs,
                                 {} if denoise else bases['atoms'], fits,
                                 atom_mapping, transfers, difference_scheme)


def calculate_directional_derivatives(measurements, plan, equilibrium=None):
    """Consume complex measured dielectric responses without reconstructing them."""
    return {atom: apply_finite_differences([epsilon for _, epsilon in values], plan.pairs[atom],
                                          plan.difference_scheme, equilibrium)
            for atom, values in measurements.items()}


def reconstruct_atomic_field(directional, plan):
    """Solve Cartesian representatives, then transfer them to equivalent atoms."""
    field = np.zeros((plan.natoms, 3, 3, 3), dtype=complex)
    for atom, slices in directional.items():
        responses = [entry[1] for entry in slices]
        if atom in plan.fits:
            tensor = apply_allowed_fit(responses, plan.fits[atom])
        else:
            basis = plan.bases[str(atom)]
            basis_responses = []
            for entry in basis['entries']:
                response = responses[entry['seed']]
                for rotation in entry['rotations']:
                    response = rotate_slice(np.asarray(rotation), response)
                basis_responses.append(response)
            fractional = np.einsum('jab,ji->abi', basis_responses, basis['inverse'])
            tensor = np.einsum('ic,abc->abi', plan.inverse_lattice, fractional)
        field[atom - 1] = tensor
    for atom, representative in plan.atom_mapping.items():
        if atom not in directional:
            field[atom - 1] = rotate_tensor3(plan.transfers[(representative, atom)], field[representative - 1])
    return field


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
    directions, responses = zip(*measured)
    pairs = prepare_pairs(directions, site_symmetry_frac, lattice, difference_scheme, atol,
                          equilibrium_available=equilibrium is not None)
    return apply_finite_differences(responses, pairs, difference_scheme, equilibrium)


def reconstruct_atomic_derivatives(measurements, lattice, site_symmetries, atom_mapping,
                                   mapping_matrices_fractional, bases=None,
                                   difference_scheme="central", equilibrium=None, denoise=False):
    """File-oriented adapter for analysis scripts; the solver uses complex arrays.

    Measurements use the XML reader's (direction, real, imaginary) records.
    The returned real/imaginary dictionaries are for legacy atomic text exports.
    Both candidate bases use the same prepared reconstruction implementation.
    """
    inverse = np.linalg.inv(lattice)
    displacements = [SimpleNamespace(atom_index=atom, fractional_vector=direction @ inverse)
                     for atom, values in measurements.items() for direction, _, _ in values]
    transfers = {pair: frac_to_cart_rotation(rotation, lattice)
                 for pair, rotation in mapping_matrices_fractional.items()}
    plan = prepare_reconstruction(displacements, lattice, site_symmetries, atom_mapping,
                                  transfers, difference_scheme, denoise, bases)
    complex_measurements = {atom: [(direction, real + 1j * imag) for direction, real, imag in values]
                            for atom, values in measurements.items()}
    reference = None if equilibrium is None else equilibrium[0] + 1j * equilibrium[1]
    directional = calculate_directional_derivatives(complex_measurements, plan, reference)
    field = reconstruct_atomic_field(directional, plan)
    return ({atom + 1: tensor.real for atom, tensor in enumerate(field)},
            {atom + 1: tensor.imag for atom, tensor in enumerate(field)})
