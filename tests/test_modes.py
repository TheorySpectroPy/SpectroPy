from pathlib import Path
import shutil
import numpy as np
import pytest
from spectropy_modes import generate, reconstruct, orbit_plan, allowed_tensor_dimension
from spectropy_structure import read_structure
from spectropy_phonons import read_gamma_modes, align_mode_atoms
from spectropy_checks import symmetry_operations, project_derivatives, mode_tensors

FIXTURE = Path(__file__).parent/"data/mos2_modes"


def test_odd_seed_and_inactive_inversion():
    angle=2*np.pi/3
    rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    mirror=np.diag([1.,-1.])
    gamma=np.array([np.eye(2),rotation,rotation@rotation,mirror,rotation@mirror,rotation@rotation@mirror])
    plan=orbit_plan(gamma)
    assert len(plan)==1 and plan[0]["minus_operation"] is not None
    assert allowed_tensor_dimension(np.array([[[1.]], [[-1.]]]),[np.eye(3),-np.eye(3)])==0


@pytest.mark.parametrize("reduced", [True,False])
def test_complex_mode_reconstruction(tmp_path,monkeypatch,reduced):
    monkeypatch.chdir(tmp_path)
    for name in ("CONTCAR","band.yaml","irreps.yaml","symmetry"):
        shutil.copy(FIXTURE/name,name)
    Path("input").write_text("1 0 0\n1 0 0\nz\nlaser_energies: 1.96\n")
    modes=None if reduced else [4,5,6,7,8]
    manifest=generate(q=.01,modes=modes,reduce=reduced)
    assert len(list(Path("mode_displacements").glob("ra_POSCAR_*"))) == (4 if reduced else 10)
    structure=read_structure()
    phonons,*_=align_mode_atoms(read_gamma_modes(complex_eigenvectors=True),structure)
    operations=symmetry_operations(structure)
    rng=np.random.default_rng(32)
    raw=rng.normal(size=(3,3,3,3))+1j*rng.normal(size=(3,3,3,3))
    raw=(raw+raw.transpose(0,1,3,2))/2
    derivative=project_derivatives(raw,operations)
    expected=mode_tensors(derivative,phonons.eigendisplacements)*abs(np.linalg.det(structure.lattice))/(4*np.pi)
    for block in manifest["groups"]:
        for seed in block["seeds"]:
            for sign,name in [(1,seed["positive"]),(-1,seed["negative"])]:
                if name is None:
                    continue
                epsilon=2*np.eye(3)+sign*np.einsum("naij,na->ij",derivative,seed["cartesian_displacement"])
                folder=Path("mode_displacements")/name
                fields = '<field>energy</field>' + ''.join(
                    f'<field>{a}{b}</field>' for a in 'xyz' for b in 'xyz')
                blocks = []
                for label, values in [('real', epsilon.real), ('imag', epsilon.imag)]:
                    rows = ''.join('<r>' + str(e) + ' ' + ' '.join(
                        f'{v:.16e}' for v in values.flat) + '</r>' for e in [0., 1.96])
                    blocks.append(f'<{label}><array>{fields}<set>{rows}</set></array></{label}>')
                (folder/'vasprun.xml').write_text(
                    '<modeling><calculation><dielectricfunction comment="density-density">'
                    + ''.join(blocks) + '</dielectricfunction></calculation></modeling>')
    reconstruct()
    result=np.load("mode_results/mode_tensors_1.96eV.npz")
    assert result["modes"].tolist()==[4,5,6,7,8]
    assert np.allclose(result["tensors"],expected[result["modes"]-1],atol=2e-6,rtol=2e-6)
    with pytest.raises(FileExistsError):
        generate(q=.01)


def test_selected_partner_expands_group(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name in ("CONTCAR","band.yaml","irreps.yaml","symmetry"):
        shutil.copy(FIXTURE/name,name)
    result=generate(q=.01,modes=[6])
    assert result["groups"][0]["modes"]==[6,7]
    assert len(list(Path("mode_displacements").glob("ra_POSCAR_*")))==1


def test_reject_nonpositive_q():
    for q in [0.,-1.,float("nan")]:
        with pytest.raises(ValueError):
            generate(q=q)


def test_mode_projector_complex_partner_noise():
    from spectropy_modes import project_mode_tensors, mode_covariance_residuals
    angle = 2*np.pi/3
    c = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    mirror = np.diag([1., -1.])
    gamma = np.array([np.eye(2), c, c@c, mirror, c@mirror, c@c@mirror])
    rotations = np.array([np.block([[g, np.zeros((2,1))], [np.zeros((1,2)), np.ones((1,1))]]) for g in gamma])
    rng = np.random.default_rng(321)
    raw = rng.normal(size=(2,3,3))+1j*rng.normal(size=(2,3,3))
    clean = project_mode_tensors(raw, gamma, rotations)
    assert np.max(mode_covariance_residuals(clean, gamma, rotations)) < 1e-12
    assert np.allclose(project_mode_tensors(clean, gamma, rotations), clean)
    assert abs(np.vdot(clean, raw-clean)) < 1e-12
    assert np.linalg.norm(raw-clean) > .1
    assert np.linalg.norm(clean.imag) > .1
    # Independent single-partner symmetrization would destroy allowed E response.
    assert np.linalg.norm(clean) > .1


def test_allowed_error_survives_zero_neumann_residual():
    """Independent analytic inversion example: symmetry is not accuracy."""
    from spectropy_modes import project_mode_tensors, mode_covariance_residuals
    rotations = np.array([np.eye(3), -np.eye(3)])
    gamma = np.ones((2, 1, 1))  # even one-dimensional mode
    exact = np.diag([1., 2., 3.])[None].astype(complex)
    allowed_error = np.diag([.2, -.4, .1])[None]*(1+2j)
    estimate = exact+allowed_error
    assert np.allclose(project_mode_tensors(estimate, gamma, rotations), estimate)
    assert np.max(mode_covariance_residuals(estimate, gamma, rotations)) == 0
    assert np.linalg.norm(estimate-exact) > .5
    # Even a measurable polarized intensity changes despite exact covariance.
    assert not np.isclose(abs(estimate[0,0,0])**2, abs(exact[0,0,0])**2)
