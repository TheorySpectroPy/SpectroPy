"""Explicit Raman-active normal-mode displacements and central-difference tensors.

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
from spectropy_symmetry import symmetry_operations, mode_representations
from spectropy_config import read_settings
from spectropy_paths import link_vasp_inputs


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


def _reference(reference):
    if reference:
        return reference
    return "ref_poscar.vasp" if Path("ref_poscar.vasp").exists() else "CONTCAR"


def generate(q=None, modes=None, reference=None,
             output=None, minimum_frequency=None):
    """Enumerate +/-Q for every Raman-active partner by default.

    Every selected mode has its own explicit positive and negative calculation.
    Explicit mode IDs bypass automatic activity/frequency filters.
    """
    settings = read_settings()
    q = settings.mode_q if q is None else q
    modes = settings.mode_indices if modes is None else modes
    output = settings.mode_displacements_directory if output is None else output
    minimum_frequency = settings.mode_minimum_frequency if minimum_frequency is None else minimum_frequency
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
    symmetry_validated = modes is None
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
    directories = [root / f'ra_POSCAR_mode{int(mode)+1}{suffix}'
                   for group, _ in planned for mode in group for suffix in ('a', 'b')]
    link_vasp_inputs(directories)
    manifest = dict(schema=2, calculator="vasp",
                    q=q, q_units="Angstrom sqrt(amu)", reference=str(Path(reference).resolve()),
                    reference_sha256=hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
                    band_sha256=hashlib.sha256(Path("band.yaml").read_bytes()).hexdigest(),
                    lattice=structure.lattice.tolist(), fractional_positions=structure.fractional_positions.tolist(),
                    symbols=structure.atom_symbols,
                    frequencies_thz=phonons.frequencies_thz.tolist(),
                    phonon_to_reference_1based=(mapping+1).tolist(),
                    origin_shift=shift.tolist(), alignment_error_angstrom=error,
                    groups=[], skipped=skipped,
                    selection="raman_active" if modes is None else "explicit",
                    minimum_frequency_cm1=minimum_frequency,
                    symmetry_validated=symmetry_validated)
    index_lines = []
    for group, gamma in planned:
        frequency = phonons.frequencies_thz[group]*33.356407889
        block = dict(modes=(group+1).tolist(),
                     ir_label=labels.get(int(group[0]+1), "---"), seeds=[])
        for number, mode in enumerate(group):
            mode_id = int(mode+1)
            stem = f"mode{mode_id}"
            displacement = q*vectors[mode]
            record = dict(label=stem, positive=f"ra_POSCAR_{stem}a",
                          negative=f"ra_POSCAR_{stem}b",
                          cartesian_displacement=displacement.tolist(),
                          mode_index=mode_id,
                          max_atomic_displacement_angstrom=float(np.max(np.linalg.norm(displacement, axis=1))))
            for sign, suffix in [(1,"a"),(-1,"b")]:
                directory = root/f"ra_POSCAR_{stem}{suffix}"
                directory.mkdir(exist_ok=True)
                fractional = structure.fractional_positions + sign*displacement @ np.linalg.inv(structure.lattice)
                write_structure(directory/"POSCAR", structure, fractional)
                # Fortran-compatible interchange files: Cartesian positions,
                # positive Q is a, negative Q is b.
                lines = [structure.comment,"1.0"]
                lines += [" ".join(f"{x:.16f}" for x in row) for row in structure.lattice]
                lines += [" ".join(structure.symbols), " ".join(map(str,structure.counts)), "cart "]
                lines += [" ".join(f"{x:.16f}" for x in row) for row in fractional @ structure.lattice]
                (root/f"POSCAR_{stem}{suffix}").write_text("\n".join(lines)+"\n")
            seed_frequency = frequency[number]
            index_lines.append(f"{mode_id:3d}{seed_frequency:10.3f}{q:10.4f}\n")
            block["seeds"].append(record)
        manifest["groups"].append(block)
    (root/"mode_index").write_text("".join(index_lines))
    (root/"mode_manifest.json").write_text(json.dumps(manifest,indent=2))
    (root/"README.txt").write_text(
        "Q units: Angstrom sqrt(amu). Physical displacement = Q * e_k/sqrt(m).\n"
        "Both signs of every selected mode are independently calculated.\n"
        "a = positive, b = negative. mode_index columns are mode, cm^-1, Q.\n"
        "Never relax a mode-displaced optical structure. This is central difference.\n")
    print(f"Generated {len(list(root.glob('ra_POSCAR_mode*')))} signed jobs for "
          f"{sum(len(g['modes']) for g in manifest['groups'])} mode partners in {root}")
    return manifest


def reconstruct(directory=None, source=None, energies=None, output=None):
    from calculate_dielectric_derivatives import read_dielectric
    settings = read_settings()
    directory = settings.mode_displacements_directory if directory is None else directory
    source = settings.dielectric_source if source is None else source
    output = settings.mode_results_directory if output is None else output
    root = Path(directory)
    data = json.loads((root/"mode_manifest.json").read_text())
    if data.get("schema") != 2 or data.get("calculator") != "vasp":
        raise ValueError("Expected a schema-2 explicit VASP mode manifest; regenerate the displacement plan")
    if not np.isfinite(data["q"]) or data["q"] <= 0:
        raise ValueError("Manifest Q must be finite and positive")
    if settings.difference_scheme != "central":
        raise ValueError("Collective mode tensors currently require difference_scheme: central")
    if settings.polarization != "specific":
        raise ValueError("Collective intensity output currently requires polarization: specific")
    volume = abs(np.linalg.det(data["lattice"]))
    mode_ids = [seed["mode_index"] for block in data["groups"] for seed in block["seeds"]]
    if not mode_ids or len(set(mode_ids)) != len(mode_ids):
        raise ValueError("Manifest must contain distinct explicit modes")
    for block in data["groups"]:
        if [s["mode_index"] for s in block["seeds"]] != block["modes"]:
            raise ValueError("Each mode must have exactly one explicit displacement pair")
        if any(not s["positive"] or not s["negative"] or s["positive"] == s["negative"] for s in block["seeds"]):
            raise ValueError("Each mode requires distinct positive and negative calculations")
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
            for seed in block["seeds"]:
                positive = dielectric(seed["positive"], energy)
                negative = dielectric(seed["negative"], energy)
                tensors[seed["mode_index"]] = volume/(4*np.pi)*(positive-negative)/(2*data["q"])
        np.savez(destination/f"raw_mode_tensors_{energy:.2f}eV.npz",
                 modes=sorted(tensors), tensors=np.asarray([tensors[m] for m in sorted(tensors)]))
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
