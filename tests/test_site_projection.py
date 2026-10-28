import tempfile
import unittest
from pathlib import Path
import numpy as np
from spectropy_config import read_settings, write_settings, Settings
from reconstruct_dielectric_derivatives import project_site_tensor, project_representatives, rotate_tensor3
from reconstruct_dielectric_derivatives import reconstruct_fractional_derivatives
from spectropy_derivatives import write_derivative_outputs
from spectropy_structure import Structure
from calculate_spectrum import read_dielectric_derivatives
import csv


class SiteProjectionTests(unittest.TestCase):
    def test_fractional_solver_projects_before_transfer(self):
        h=np.diag([-1.,-1,1])
        g=np.diag([-1.,1,-1])
        raw=np.random.default_rng(5).normal(size=(3,3,3))
        measurements={1:[]}
        for axis in np.eye(3):
            for sign in (1,-1):
                d=.01*sign*axis
                eps=np.eye(3)+np.einsum('abi,i->ab',raw,d)
                measurements[1].append((d,eps,np.zeros((3,3))))
        sites={1:{'rotations':np.array([np.eye(3),h])}}
        mapping={1:1,2:1}
        kwargs=dict(measurements=measurements,lattice=np.eye(3),site_symmetries=sites,
                    atom_mapping=mapping,mapping_matrices_fractional={(1,2):g})
        uncorrected,_=reconstruct_fractional_derivatives(**kwargs)
        corrected,_=reconstruct_fractional_derivatives(**kwargs,site_symmetry_projection=True)
        np.testing.assert_allclose(uncorrected[1],raw,atol=1e-13)
        expected=project_site_tensor(raw,[np.eye(3),h])
        np.testing.assert_allclose(corrected[1],expected,atol=1e-13)
        np.testing.assert_allclose(corrected[2],rotate_tensor3(g,expected),atol=1e-13)

    def test_full_relations_complex_idempotence(self):
        r = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        group = [np.linalg.matrix_power(r, n) for n in range(4)]
        rng = np.random.default_rng(2)
        raw = rng.normal(size=(3,3,3)) + 1j*rng.normal(size=(3,3,3))
        out = project_site_tensor(raw, group)
        np.testing.assert_allclose(out, project_site_tensor(out, group), atol=1e-14)
        for g in group:
            np.testing.assert_allclose(out, rotate_tensor3(g,out), atol=1e-14)
        self.assertAlmostEqual(abs(out[0,0,2]-out[1,1,2]),0.)
        self.assertLessEqual(np.linalg.norm(out),np.linalg.norm(raw))
        self.assertAlmostEqual(abs(np.vdot(out,raw-out)),0.)

    def test_route_independence_and_field_equivalence(self):
        h = np.diag([-1.,-1,1])
        g = np.diag([-1.,1,-1])
        d = np.random.default_rng(3).normal(size=(3,3,3))
        p = project_site_tensor(d,[np.eye(3),h])
        np.testing.assert_allclose(rotate_tensor3(g,p),rotate_tensor3(g@h,p),atol=1e-14)
        # Two-site orbit: h fixes each site; g swaps sites. Start from copies
        # transferred from the raw representative and compare the group average.
        field = [d,rotate_tensor3(g,d)]
        expected = [p,rotate_tensor3(g,p)]
        for a in range(2):
            projected = (field[a]+rotate_tensor3(h,field[a])+rotate_tensor3(g,field[1-a])+rotate_tensor3(g@h,field[1-a]))/4
            np.testing.assert_allclose(projected,expected[a],atol=1e-14)

    def test_reject_extra_measured_atoms(self):
        d = np.zeros((3,3,3))
        with self.assertRaises(ValueError):
            project_representatives({1:d,2:d},{1:d,2:d},{},np.eye(3),{1:1,2:1})

    def test_identity_and_invalid_group(self):
        d=np.ones((3,3,3))
        np.testing.assert_array_equal(project_site_tensor(d,[np.eye(3)]),d)
        with self.assertRaises(ValueError):
            project_site_tensor(d,[2*np.eye(3)])
        with self.assertRaises(ValueError):
            project_site_tensor(d,[])

    def test_setting_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'input'
            self.assertFalse(read_settings(p).site_symmetry_projection)
            settings=Settings(incident_polarization=(1,0,0),scattered_polarization=(1,0,0),surface_normal='z',site_symmetry_projection=True)
            write_settings(p,settings)
            self.assertEqual(read_settings(p),settings)
            p.write_text('site_symmetry_projection: maybe\n')
            with self.assertRaises(ValueError):
                read_settings(p)

    def test_raw_corrected_outputs_and_residuals(self):
        structure=Structure('test',np.diag([2.,2.,3.]),('B','N'),np.array([1,1]),
                            np.array([[0.,0.,0.],[.5,.5,.5]]))
        raw_r={1:np.ones((3,3,3)),2:2*np.ones((3,3,3))}
        raw_i={1:np.ones((3,3,3)),2:np.zeros((3,3,3))}
        corrected_r={1:np.zeros((3,3,3)),2:raw_r[2].copy()}
        corrected_i={1:np.zeros((3,3,3)),2:raw_i[2].copy()}
        with tempfile.TemporaryDirectory() as tmp:
            rows=write_derivative_outputs(0.,structure,(raw_r,raw_i),(corrected_r,corrected_i),[1],tmp)
            self.assertAlmostEqual(rows[0]['residual'],np.sqrt(54))
            self.assertEqual(rows[1]['residual'],0.)
            self.assertTrue(rows[0]['representative'])
            self.assertFalse(rows[1]['representative'])
            expected=np.sqrt(54)*12/(4*np.pi)
            self.assertAlmostEqual(rows[0]['stored_units_residual'],expected)
            a,b=read_dielectric_derivatives(Path(tmp)/'raw_epsilon_derivative_0.00',2)
            self.assertAlmostEqual(a[0,0,0,0],12/(4*np.pi))
            self.assertAlmostEqual(b[0,0,0,0],12/(4*np.pi))
            c,_=read_dielectric_derivatives(Path(tmp)/'epsilon_derivative_0.00',2)
            self.assertAlmostEqual(c[0,0,0,0],0.)
            with (Path(tmp)/'site_projection_residuals_0.00.csv').open() as file:
                self.assertEqual(len(list(csv.DictReader(file))),2)

if __name__ == '__main__':
    unittest.main()
