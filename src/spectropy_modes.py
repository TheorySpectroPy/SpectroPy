"""Direct Raman-active normal-mode displacements; optional experimental reduction.

Q is in Angstrom sqrt(amu), as in Liangbo's raman_dis_mode.f90.
No Fortran utility or atomic-derivative field is needed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import yaml

from spectropy_structure import read_structure, write_structure
from spectropy_phonons import read_gamma_modes, align_mode_atoms
from spectropy_checks import symmetry_operations, mode_representations


def mode_covariance_residuals(tensors, representations, rotations):
    """Absolute Neumann RMS residual for each partner, including complex parts."""
    values = np.asarray(tensors)
    gamma = np.asarray(representations)
    rotations = np.asarray(rotations)
    if values.ndim != 3 or values.shape[1:] != (3, 3) or gamma.shape != (len(rotations), len(values), len(values)):
        raise ValueError("Incompatible mode tensors and partner representations")
    if not len(rotations) or not all(np.isfinite(x).all() for x in (values, gamma, rotations)):
        raise ValueError("Finite tensors and a nonempty symmetry group required")
    differences = np.asarray([
        np.einsum("ia,pab,jb->pij", rotation, values, rotation)
        - np.einsum("pk,pij->kij", matrix, values)
        for matrix, rotation in zip(gamma, rotations)])
    return np.sqrt(np.mean(np.sum(abs(differences)**2, axis=(2, 3)), axis=0))


def project_mode_tensors(tensors, representations, rotations):
    """Orthogonal group projector on a COMPLETE real orthonormal partner basis.

    The caller must supply the complete group, with matching representation
    matrices (validated from eigenvectors), not a list of generators.
    This minimizes the total squared covariance mismatch, not a sum of norms.
    """
    values = np.asarray(tensors)
    gamma = np.real_if_close(np.asarray(representations))
    rotations = np.asarray(rotations)
    mode_covariance_residuals(values, gamma, rotations)  # shape/finite validation
    if np.iscomplexobj(gamma) or not np.allclose(gamma @ gamma.transpose(0, 2, 1), np.eye(len(values)), atol=1e-7):
        raise ValueError("A real orthonormal partner representation is required")
    if not np.allclose(rotations @ rotations.transpose(0, 2, 1), np.eye(3), atol=1e-7):
        raise ValueError("Cartesian symmetry rotations must be orthogonal")
    corrected = sum(np.einsum("kp,ia,pab,jb->kij", matrix, rotation, values, rotation)
                    for matrix, rotation in zip(gamma, rotations))/len(rotations)
    if np.max(mode_covariance_residuals(corrected, gamma, rotations)) > 1e-7*max(1., np.linalg.norm(values)):
        raise ValueError("Supplied operations do not yield a covariant projector")
    return corrected


def allowed_tensor_dimension(representations, rotations, tolerance=1e-7):
    """Dimension of symmetric rank-two intertwiners, not a label lookup."""
    d = representations.shape[1]
    constraints = []
    for gamma, rotation in zip(representations, rotations):
        left = np.kron(np.eye(d), np.kron(rotation, rotation))
        right = np.kron(gamma.T, np.eye(9))
        constraints.append(left-right)
    symmetry = np.zeros((9,9))
    for a in range(3):
        for b in range(3):
            symmetry[3*a+b,3*a+b] += .5
            symmetry[3*a+b,3*b+a] += .5
    constraints.append(np.kron(np.eye(d), np.eye(9)-symmetry))
    singular = np.linalg.svd(np.vstack(constraints), compute_uv=False)
    return int(9*d - np.count_nonzero(singular > tolerance))


def orbit_plan(representations, reuse_sign=True, tolerance=1e-7):
    """Greedy full-rank orbit coverage; no global-minimum claim."""
    gamma = np.real_if_close(representations)
    if np.iscomplexobj(gamma):
        raise ValueError("Complex partner bases require realification; not supported")
    d = gamma.shape[1]
    candidates = [np.eye(d)[:,i] for i in range(d)]
    if reuse_sign:
        for matrix in gamma:
            # SVD nullspace of Gamma+I yields real -1 eigenvectors even for
            # nonsymmetric representation matrices.
            _, singular, vh = np.linalg.svd(matrix+np.eye(d))
            for index in np.flatnonzero(singular < tolerance):
                candidates.append(vh[index])
    selected, basis = [], np.empty((d,0))
    while np.linalg.matrix_rank(basis, tol=tolerance) < d:
        best = None
        old_rank = np.linalg.matrix_rank(basis, tol=tolerance)
        for vector in candidates:
            vector = np.asarray(vector, dtype=float)
            vector /= np.linalg.norm(vector)
            if vector[np.argmax(abs(vector))] < 0:
                vector = -vector
            orbit = np.einsum("gij,j->gi", gamma, vector).T
            minus = next((i for i,v in enumerate(orbit.T)
                          if np.linalg.norm(v+vector) < tolerance), None) if reuse_sign else None
            cost = 1 if minus is not None else 2
            gain = np.linalg.matrix_rank(np.column_stack([basis,orbit]), tol=tolerance)-old_rank
            score = (gain/cost, gain, -cost)
            if gain and (best is None or score > best[0]):
                best = score, vector.copy(), minus, orbit
        if best is None:
            raise ValueError("Could not span partner space with symmetry orbits")
        _, vector, minus, orbit = best
        selected.append(dict(coefficients=vector.tolist(), minus_operation=minus))
        basis = np.column_stack([basis,orbit])
    return selected


def _reference(reference):
    if reference:
        return reference
    return "ref_poscar.vasp" if Path("ref_poscar.vasp").exists() else "CONTCAR"


def generate(q=.01, modes=None, reference=None, reduce=False,
             reuse_sign=True, output="mode_displacements", minimum_frequency=1.0):
    """Enumerate +/-Q for every Raman-active partner by default.

    Symmetry selects active groups but does not reconstruct measurements unless
    reduce=True. Explicit mode IDs bypass automatic activity/frequency filters.
    """
    if not np.isfinite(q) or q <= 0:
        raise ValueError("Q must be finite and positive (Angstrom sqrt(amu))")
    if not np.isfinite(minimum_frequency) or minimum_frequency < 0:
        raise ValueError("Minimum frequency must be finite and nonnegative (cm^-1)")
    reference = _reference(reference)
    structure = read_structure(reference)
    phonons = read_gamma_modes(complex_eigenvectors=True)
    phonons, mapping, shift, error = align_mode_atoms(phonons, structure)
    if not np.isfinite(phonons.eigendisplacements).all() or not np.isfinite(phonons.frequencies_thz).all() or np.any(phonons.masses <= 0):
        raise ValueError("Finite Gamma eigenvectors/frequencies and positive masses required")
    vectors = np.real_if_close(phonons.eigendisplacements)
    if np.iscomplexobj(vectors):
        if np.max(abs(vectors.imag)) > 1e-8:
            raise ValueError("Mode displacements require a real Gamma eigenvector basis")
        vectors = vectors.real
    nmodes = len(vectors)
    requested = set(range(nmodes)) if modes is None else {int(m)-1 for m in modes}
    if not requested or min(requested) < 0 or max(requested) >= nmodes:
        raise ValueError("Requested mode index outside 1..3N")
    # Activity selection requires the actual crystal representation even when
    # every partner and both signs will be computed independently.
    symmetry_validated = reduce or modes is None
    labels = {}
    if symmetry_validated:
        entries = yaml.safe_load(Path("irreps.yaml").read_text())
        if not np.allclose(entries["q-position"], 0, atol=1e-8):
            raise ValueError("irreps.yaml must describe Gamma")
        groups = [np.asarray(e["band_indices"],dtype=int)-1 for e in entries["normal_modes"]]
        labels = {int(i): e.get("ir_label", "---")
                  for e in entries["normal_modes"] for i in e["band_indices"]}
        operations = symmetry_operations(structure)
        representations, leakages = mode_representations(phonons, operations, groups, tolerance=1e-5)
        rotations = [r for _,r in operations]
    else:
        groups = [np.array([m]) for m in range(nmodes)]
        representations = [np.ones((1,1,1)) for _ in groups]
        rotations = [np.eye(3)]
    planned, skipped = [], []
    for group, gamma in zip(groups, representations):
        if not requested.intersection(group.tolist()):
            continue
        frequency = phonons.frequencies_thz[group]*33.356407889
        if np.ptp(frequency) > 1e-2:
            raise ValueError("Partner group not frequency-degenerate to 0.01 cm^-1")
        if modes is None and (np.min(frequency) <= minimum_frequency or
                              allowed_tensor_dimension(gamma, rotations) == 0):
            skipped.extend((group+1).tolist())
            continue
        planned.append((group, gamma))
    if not planned:
        raise ValueError("No Raman-active modes above the frequency cutoff; no displacements generated")
    root = Path(output)
    if root.exists():
        raise FileExistsError(f"Use a fresh mode-displacement directory: {root}")
    root.mkdir(parents=True)
    manifest = dict(schema=1, calculator="vasp",
                    q=q, q_units="Angstrom sqrt(amu)", reference=str(Path(reference).resolve()),
                    reference_sha256=hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
                    band_sha256=hashlib.sha256(Path("band.yaml").read_bytes()).hexdigest(),
                    lattice=structure.lattice.tolist(), fractional_positions=structure.fractional_positions.tolist(),
                    symbols=structure.atom_symbols, rotations=[r.tolist() for r in rotations],
                    frequencies_thz=phonons.frequencies_thz.tolist(),
                    phonon_to_reference_1based=(mapping+1).tolist(),
                    origin_shift=shift.tolist(), alignment_error_angstrom=error,
                    groups=[], skipped=skipped, reduction=reduce,
                    selection="raman_active" if modes is None else "explicit",
                    minimum_frequency_cm1=minimum_frequency,
                    symmetry_validated=symmetry_validated)
    index_lines = []
    for group, gamma in planned:
        frequency = phonons.frequencies_thz[group]*33.356407889
        seeds = (orbit_plan(gamma, reuse_sign) if reduce else
                 [dict(coefficients=row.tolist(), minus_operation=None) for row in np.eye(len(group))])
        block = dict(modes=(group+1).tolist(), representations=np.asarray(gamma).real.tolist(),
                     ir_label=labels.get(int(group[0]+1), "---"), seeds=[])
        for number, seed in enumerate(seeds):
            # Mode index plus optional seed suffix; JSON is authoritative for
            # rotated partner combinations, not the filename alone.
            mode_id = int(group[0]+1 if reduce else group[number]+1)
            stem = f"mode{mode_id}" + (f"_seed{number+1}" if reduce and len(seeds)>1 else "")
            displacement = q*np.einsum("p,pni->ni", seed["coefficients"], vectors[group])
            record = dict(**seed, label=stem, positive=f"ra_POSCAR_{stem}a",
                          negative=None if seed["minus_operation"] is not None else f"ra_POSCAR_{stem}b",
                          cartesian_displacement=displacement.tolist(),
                          mode_index=mode_id,
                          max_atomic_displacement_angstrom=float(np.max(np.linalg.norm(displacement, axis=1))))
            for sign, suffix in [(1,"a"),(-1,"b")]:
                if sign == -1 and record["negative"] is None:
                    continue
                directory = root/f"ra_POSCAR_{stem}{suffix}"
                directory.mkdir()
                fractional = structure.fractional_positions + sign*displacement @ np.linalg.inv(structure.lattice)
                write_structure(directory/"POSCAR", structure, fractional)
                # Fortran-compatible interchange files: Cartesian positions,
                # positive Q is a, negative Q is b.
                lines = [structure.comment,"1.0"]
                lines += [" ".join(f"{x:.16f}" for x in row) for row in structure.lattice]
                lines += [" ".join(structure.symbols), " ".join(map(str,structure.counts)), "cart "]
                lines += [" ".join(f"{x:.16f}" for x in row) for row in fractional @ structure.lattice]
                (root/f"POSCAR_{stem}{suffix}").write_text("\n".join(lines)+"\n")
            seed_frequency = np.mean(frequency) if reduce else frequency[number]
            index_lines.append(f"{mode_id:3d}{seed_frequency:10.3f}{q:10.4f}\n")
            block["seeds"].append(record)
        manifest["groups"].append(block)
    (root/"mode_index").write_text("".join(index_lines))
    (root/"mode_manifest.json").write_text(json.dumps(manifest,indent=2))
    (root/"README.txt").write_text(
        "Q units: Angstrom sqrt(amu). Physical displacement = Q * sum(c_p e_p/sqrt(m)).\n"
        f"Partner/sign reduction: {reduce}. Without reduction each seed is one original mode,\n"
        "and both signs of every selected mode are independently calculated.\n"
        "a = positive, b = negative. mode_index has representative mode, cm^-1, Q,\n"
        "matching Liangbo's columns. In reduced runs a seed can be a linear\n"
        "combination of degenerate partners: read mode_manifest.json coefficients;\n"
        "the filename alone is NOT a specification of its eigenvector.\n"
        "Omitted negative files are supplied by the recorded exact symmetry operation.\n"
        "Never relax a mode-displaced optical structure. This is central difference.\n")
    (root/"runRA").write_text(
        "#!/bin/bash\nset -euo pipefail\n# Directories already exist. Stage only calculator inputs, never run DFT here.\n"
        'cd "$(dirname "$0")"\n'
        "for d in ra_POSCAR_mode*; do\n"
        '  [ -d "$d" ] || continue\n'
        '  for f in INCAR KPOINTS POTCAR; do\n'
        '    [ ! -f "../$f" ] || cp "../$f" "$d/"\n'
        "  done\n"
        "done\n")
    print(f"Generated {len(list(root.glob('ra_POSCAR_mode*')))} signed jobs for "
          f"{sum(len(g['modes']) for g in manifest['groups'])} mode partners in {root}")
    return manifest


def reconstruct(directory="mode_displacements", source="vasprun.xml", energies=None, output="mode_results", denoise=False):
    from calculate_dielectric_derivatives import read_dielectric
    from spectropy_config import read_settings
    root = Path(directory)
    data = json.loads((root/"mode_manifest.json").read_text())
    if data["calculator"] != "vasp":
        raise ValueError("Only VASP mode manifests are supported")
    settings = read_settings()
    if settings.difference_scheme != "central":
        raise ValueError("Collective mode tensors currently require difference_scheme: central")
    if settings.polarization != "specific":
        raise ValueError("Collective intensity output currently requires polarization: specific")
    rotations = np.asarray(data["rotations"])
    volume = abs(np.linalg.det(data["lattice"]))
    projection_groups = None
    if denoise:
        if data["calculator"] != "vasp":
            raise ValueError("Mode denoising is currently validated on the VASP development path only")
        if hashlib.sha256(Path("band.yaml").read_bytes()).hexdigest() != data["band_sha256"]:
            raise ValueError("Projection band.yaml differs from the displacement basis")
        structure = read_structure(data["reference"])
        phonons, *_ = align_mode_atoms(read_gamma_modes(complex_eigenvectors=True), structure)
        entries = yaml.safe_load(Path("irreps.yaml").read_text())
        if not np.allclose(entries["q-position"], 0, atol=1e-8):
            raise ValueError("Mode denoising requires Gamma irreps")
        projection_groups = [np.asarray(e["band_indices"], dtype=int)-1 for e in entries["normal_modes"]]
        operations = symmetry_operations(structure)
        projection_gamma, _ = mode_representations(phonons, operations, projection_groups, tolerance=1e-5)
        projection_rotations = np.asarray([r for _, r in operations])
    destination = Path(output)
    if destination.exists():
        raise FileExistsError(f"Use a fresh output directory: {destination}")
    destination.mkdir(parents=True)
    def dielectric(name, energy):
        folder = root/name
        actual = read_structure(folder/"POSCAR")
        seed = next(seed for block in data["groups"] for seed in block["seeds"]
                    if name in (seed["positive"],seed["negative"]))
        sign = 1 if name == seed["positive"] else -1
        expected = np.asarray(data["fractional_positions"]) + sign*np.asarray(seed["cartesian_displacement"]) @ np.linalg.inv(data["lattice"])
        delta = actual.fractional_positions-expected
        if actual.atom_symbols != data["symbols"] or not np.allclose(actual.lattice,data["lattice"],atol=1e-9,rtol=0) or np.max(np.linalg.norm((delta-np.rint(delta)) @ actual.lattice,axis=1)) > 1e-7:
            raise ValueError(f"Mode-displacement geometry/order mismatch in {folder}")
        real, imag = read_dielectric(folder/"vasprun.xml", energy, source)
        if real is None:
            raise ValueError(f"Missing dielectric data in {folder}")
        return real+1j*imag
    for energy in settings.laser_energies if energies is None else energies:
        tensors = {}
        for block in data["groups"]:
            gamma = np.asarray(block["representations"])
            directions, responses = [], []
            for seed in block["seeds"]:
                positive = dielectric(seed["positive"], energy)
                operation = seed["minus_operation"]
                negative = (rotations[operation] @ positive @ rotations[operation].T
                            if operation is not None else dielectric(seed["negative"], energy))
                derivative = (positive-negative)/(2*data["q"])
                # Without reduction, use ONLY the directly calculated +/- pair.
                # Symmetry used for activity selection must not average raw data.
                images = (zip(gamma, rotations) if data["reduction"] else
                          [(np.eye(len(block["modes"])), np.eye(3))])
                for matrix, rotation in images:
                    directions.append(matrix @ seed["coefficients"])
                    responses.append((rotation @ derivative @ rotation.T).reshape(9))
            directions = np.asarray(directions).T
            if np.linalg.matrix_rank(directions, tol=1e-7) != len(block["modes"]):
                raise ValueError("Insufficient mode-orbit rank")
            partners = np.asarray(responses).T @ np.linalg.pinv(directions)
            for index, mode in enumerate(block["modes"]):
                tensors[mode] = volume/(4*np.pi)*partners[:,index].reshape(3,3)
        raw_tensors = {m: tensor.copy() for m, tensor in tensors.items()}
        if denoise:
            import csv
            rows = []
            for group, gamma in zip(projection_groups, projection_gamma):
                ids = (group+1).tolist()
                if not set(ids).intersection(tensors):
                    continue
                if not set(ids).issubset(tensors):
                    raise ValueError("Denoising requires every member of each selected partner group")
                values = np.asarray([tensors[m] for m in ids])
                corrected = project_mode_tensors(values, gamma, projection_rotations)
                before = mode_covariance_residuals(values, gamma, projection_rotations)
                after = mode_covariance_residuals(corrected, gamma, projection_rotations)
                for m, raw, clean, rho0, rho1 in zip(ids, values, corrected, before, after):
                    tensors[m] = clean
                    rows.append([m, np.linalg.norm(raw), np.linalg.norm(raw-clean), rho0, rho1])
            with (destination/f"mode_projection_residuals_{energy:.2f}eV.csv").open("w") as stream:
                writer = csv.writer(stream)
                writer.writerow(["mode", "raw_norm", "removed_norm", "raw_neumann_absolute", "corrected_neumann_absolute"])
                writer.writerows(rows)
            np.savez(destination/f"corrected_mode_tensors_{energy:.2f}eV.npz",
                     modes=sorted(tensors), tensors=np.asarray([tensors[m] for m in sorted(tensors)]))
        np.savez(destination/f"raw_mode_tensors_{energy:.2f}eV.npz",
                 modes=sorted(raw_tensors), tensors=np.asarray([raw_tensors[m] for m in sorted(raw_tensors)]))
        np.savez(destination/f"mode_tensors_{energy:.2f}eV.npz",
                 modes=sorted(tensors), tensors=np.asarray([tensors[m] for m in sorted(tensors)]))
        import csv
        with (destination/f"mode_covariance_residuals_{energy:.2f}eV.csv").open("w") as stream:
            writer=csv.writer(stream)
            writer.writerow(["mode","frequency_cm-1","absolute_residual"])
            for block in data["groups"]:
                partners=np.asarray([tensors[m] for m in block["modes"]])
                for index,mode in enumerate(block["modes"]):
                    differences=[]
                    for gamma,rotation in zip(block["representations"],rotations):
                        differences.append(rotation @ tensors[mode] @ rotation.T -
                                           np.einsum("p,pij->ij",np.asarray(gamma)[:,index],partners))
                    rho=np.sqrt(np.mean(np.sum(abs(np.asarray(differences))**2,axis=(1,2))))
                    # An identity-only explicit-mode plan cannot measure full
                    # crystal covariance; do not misreport its trivial zero.
                    if not data.get("symmetry_validated", data["reduction"]):
                        rho = float("nan")
                    writer.writerow([mode,data["frequencies_thz"][mode-1]*33.356407889,rho])
        with (destination/f"Raman_tensor_{energy:.2f}eV").open("w") as stream, \
             (destination/f"Raman_intensity_complex_{energy:.2f}eV").open("w") as intensity:
            stream.write("# Mode Freq(THz) Freq(cm-1) Raman Tensor (Real + i*Imaginary)\n")
            for mode, tensor in sorted(tensors.items()):
                frequency = data["frequencies_thz"][mode-1]
                stream.write(f"{mode:5d} {frequency:.12e} {frequency*33.356407889:.12e}\n")
                for row in tensor:
                    stream.write(" ".join(f"{z.real:.16e}{z.imag:+.16e}j" for z in row)+"\n")
                stream.write("\n")
                if settings.incident_polarization is not None and settings.scattered_polarization is not None:
                    value = abs(np.asarray(settings.scattered_polarization) @ tensor @ settings.incident_polarization)**2
                    intensity.write(f"{frequency*33.356407889:.12e} {value:.16e} ---\n")
    print(f"Reconstructed collective mode tensors in {destination}")
