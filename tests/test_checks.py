import csv
from pathlib import Path
import tempfile
import unittest

import numpy as np
import yaml

from spectropy_checks import angular_intensities, project_derivatives, projection_residuals
from spectropy_checks import align_mode_atoms, mode_representations, mode_covariance_residuals
from spectropy_checks import mode_tensors
from spectropy_structure import Structure
from spectropy_cli import main
from spectropy_phonons import read_gamma_modes
from spectropy_phonons import GammaModes


class ModeCovarianceChecks(unittest.TestCase):
    def setUp(self):
        self.operations = [(np.array([0]), np.eye(3)),
                           (np.array([0]), np.diag([-1., -1., 1.]))]
        self.phonons = GammaModes(np.array([1., 1., 2.]), np.eye(3).reshape(3, 1, 3), np.ones(1), 1)
        self.groups = [np.array([0, 1]), np.array([2])]
        self.tensors = np.array([[[0, 0, 1], [0, 0, 0], [1, 0, 0]],
                                 [[0, 0, 0], [0, 0, 1], [0, 1, 0]], np.eye(3)], dtype=complex)

    def test_exact_covariance_and_zero_tensors(self):
        reps, _ = mode_representations(self.phonons, self.operations, self.groups)
        for tensors in (self.tensors, self.tensors * 0):
            modes, groups = mode_covariance_residuals(tensors, self.operations, self.groups, reps)
            np.testing.assert_allclose(modes["absolute_residual"], 0, atol=1e-14)
            self.assertNotIn('relative_residual', modes)
            np.testing.assert_allclose(modes['comparison_norm'], np.maximum(modes['lhs_norm'], modes['rhs_norm']))
            np.testing.assert_allclose(groups["absolute_residual"], 0, atol=1e-14)

    def test_detects_forbidden_component(self):
        reps, _ = mode_representations(self.phonons, self.operations, self.groups)
        self.tensors[0, 0, 0] = 2j
        modes, _ = mode_covariance_residuals(self.tensors, self.operations, self.groups, reps)
        self.assertAlmostEqual(modes["squared_residual"][0], 8)
        self.assertGreater(modes["absolute_residual"][0], 0)

    def test_complex_partner_basis_invariance_and_nondiagonal_gamma(self):
        # A quarter-turn mixes x and y, testing Gamma_pk orientation, not just signs.
        self.operations.append((np.array([0]), np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])))
        reps, _ = mode_representations(self.phonons, self.operations, self.groups)
        exact, _ = mode_covariance_residuals(self.tensors, self.operations, self.groups, reps)
        np.testing.assert_allclose(exact["absolute_residual"], 0, atol=1e-14)
        self.tensors[0, 0, 0] = 2j
        _, before = mode_covariance_residuals(self.tensors, self.operations, self.groups, reps)
        change = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
        mixed = self.tensors.copy()
        mixed[:2] = np.einsum("kp,pij->kij", change, self.tensors[:2])
        eig = self.phonons.eigendisplacements.astype(complex).copy()
        eig[:2] = np.einsum("kp,pna->kna", change, eig[:2])
        phonons = GammaModes(self.phonons.frequencies_thz, eig, self.phonons.masses, 1)
        new_reps, _ = mode_representations(phonons, self.operations, self.groups)
        _, after = mode_covariance_residuals(mixed, self.operations, self.groups, new_reps)
        np.testing.assert_allclose(before["absolute_residual"], after["absolute_residual"], atol=1e-14)
        np.testing.assert_allclose(before["comparison_norm"], after["comparison_norm"], atol=1e-14)

    def test_rejects_nonclosed_partner_groups(self):
        self.operations.append((np.array([0]), np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])))
        with self.assertRaisesRegex(ValueError, "not a closed"):
            mode_representations(self.phonons, self.operations, [np.array([i]) for i in range(3)])

    def test_atom_reordering_and_origin_alignment(self):
        structure = Structure("test", np.eye(3), ("B", "N"), np.array([1, 2]),
                              np.array([[0.1, 0, 0], [0.4, 0, 0], [0.8, 0, 0]]))
        order = np.array([2, 0, 1])
        shift = np.array([0.2, 0.1, 0])
        masses = np.array([14., 10., 14.])
        eig = np.eye(9).reshape(9, 3, 3) / np.sqrt(masses)[None, :, None]
        phonons = GammaModes(np.arange(9), eig, masses, 3)
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "band.yaml"
            path.write_text(yaml.safe_dump({"lattice": np.eye(3).tolist(), "points": [
                {"symbol": structure.atom_symbols[i], "coordinates": (structure.fractional_positions[i] - shift).tolist()}
                for i in order]}))
            aligned, mapping, found_shift, error = align_mode_atoms(phonons, structure, path)
        np.testing.assert_array_equal(mapping, order)
        np.testing.assert_allclose(found_shift, shift)
        np.testing.assert_allclose(aligned.eigendisplacements[:, mapping], eig)
        np.testing.assert_allclose(aligned.masses[mapping], masses)
        self.assertLess(error, 1e-14)

    def test_registered_contraction_is_invariant_to_cyclic_atom_order(self):
        # Four-cycle, unequal masses, complex eigenvectors and complex atomic
        # responses: catch permutation direction and mass-weighting regressions.
        rng = np.random.default_rng(7)
        positions = np.array([[.11, .13, .17], [.37, .19, .29], [.53, .41, .23], [.79, .61, .43]])
        structure = Structure('test', np.eye(3), ('B',), np.array([4]), positions)
        order = np.array([1, 2, 3, 0])
        shift = np.array([.21, .07, .03])
        masses = np.array([10., 11., 12., 13.])
        physical = (rng.normal(size=(12, 4, 3)) + 1j * rng.normal(size=(12, 4, 3))) / np.sqrt(masses)[None, :, None]
        source = GammaModes(np.arange(12), physical[:, order], masses[order], 4)
        derivatives = rng.normal(size=(4, 3, 3, 3)) + 1j * rng.normal(size=(4, 3, 3, 3))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'band.yaml'
            path.write_text(yaml.safe_dump(dict(lattice=np.eye(3).tolist(), points=[dict(symbol='B', coordinates=(positions[i] - shift).tolist()) for i in order])))
            aligned, _, _, _ = align_mode_atoms(source, structure, path)
        np.testing.assert_allclose(aligned.eigendisplacements, physical, atol=1e-14)
        np.testing.assert_allclose(aligned.masses, masses)
        np.testing.assert_allclose(mode_tensors(derivatives, aligned.eigendisplacements), mode_tensors(derivatives, physical), atol=1e-14)
        self.assertGreater(np.linalg.norm(mode_tensors(derivatives, source.eigendisplacements) - mode_tensors(derivatives, physical)), 1)


class TensorChecks(unittest.TestCase):
    def test_inversion_forbids_rank_three_response(self):
        tensor = np.arange(27).reshape(1, 3, 3, 3).astype(complex) * (1 + 2j)
        operations = [(np.array([0]), np.eye(3)), (np.array([0]), -np.eye(3))]
        np.testing.assert_allclose(project_derivatives(tensor, operations), 0)

    def test_equivalent_atom_permutation_and_idempotence(self):
        tensors = np.arange(54).reshape(2, 3, 3, 3).astype(float)
        operations = [(np.array([0, 1]), np.eye(3)), (np.array([1, 0]), -np.eye(3))]
        projected = project_derivatives(tensors, operations)
        np.testing.assert_allclose(projected[0], (tensors[0] - tensors[1]) / 2)
        np.testing.assert_allclose(projected[1], -projected[0])
        np.testing.assert_allclose(project_derivatives(projected, operations), projected)

    def test_residual_zero_and_symmetric_normalization(self):
        original = np.array([np.eye(3), np.zeros((3, 3))], dtype=complex)
        projected = -original
        result = projection_residuals(original, projected)
        np.testing.assert_allclose(result["relative_residual"], [2, 0])
        np.testing.assert_allclose(result["relative_residual"],
                                   projection_residuals(projected, original)["relative_residual"])
        np.testing.assert_allclose(result["squared_residual"], [12, 0])

    def test_a1_parallel_circle_crossed_zero(self):
        tensor = np.diag([1 + 2j, 1 + 2j, 3 - 4j])
        angles = np.arange(361)
        np.testing.assert_allclose(angular_intensities(tensor, angles), 5, atol=1e-14)
        np.testing.assert_allclose(angular_intensities(tensor, angles, scan="crossed"), 0, atol=1e-28)

    def test_complex_interference_and_global_phase(self):
        tensor = np.diag([1, 1j, 0])
        values = angular_intensities(tensor, [0, 45, 90, 180])
        np.testing.assert_allclose(values, [[1, 0.5, 1, 1]])
        np.testing.assert_allclose(values, angular_intensities(tensor * np.exp(0.7j), [0, 45, 90, 180]))

    def test_degenerate_incoherent_sum_is_basis_independent(self):
        tensors = np.array([np.diag([1, -1, 0]), [[0, 1, 0], [1, 0, 0], [0, 0, 0]]])
        angles = np.arange(361)
        summed = angular_intensities(tensors, angles).sum(axis=0)
        np.testing.assert_allclose(summed, 1, atol=1e-14)
        change = np.array([[np.cos(0.3), -np.sin(0.3)], [np.sin(0.3), np.cos(0.3)]])
        mixed = np.einsum("mn,nij->mij", change, tensors)
        np.testing.assert_allclose(angular_intensities(mixed, angles).sum(axis=0), summed)

    def test_invalid_diagnostic_inputs(self):
        with self.assertRaises(ValueError):
            projection_residuals(np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), 0)
        with self.assertRaises(ValueError):
            angular_intensities(np.eye(3), [float("nan")])


class CommandChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)
        modes = []
        for index in range(3):
            pairs = np.zeros((1, 3, 2))
            pairs[0, index, 0] = 1
            modes.append({"frequency": float(index + 1), "eigenvector": pairs.tolist()})
        data = {"natom": 1, "lattice": np.eye(3).tolist(),
                "points": [{"mass": 1, "symbol": "B", "coordinates": [0, 0, 0]}],
                "phonon": [{"q-position": [0, 0, 0], "band": modes}]}
        (self.path / "band.yaml").write_text(yaml.safe_dump(data))
        (self.path / "input").write_text("1 0 0\n1 0 0\nz\nlaser_energies: 0.00\n")
        (self.path / "epsilon_derivative_0.00").write_text(
            "! Real Part\nAtom 1\n"
            "1 0 0 0 0 0 0 0 0\n0 1 0 0 0 0 0 0 0\n0 0 1 0 0 0 0 0 0\n"
            "! The Imaginary Part\nAtom 1\n"
            "0 0 0 0 0 0 0 0 0\n0 0 0 0 0 0 0 0 0\n0 0 0 0 0 0 0 0 0\n"
        )

    def tearDown(self):
        self.temp.cleanup()

    def test_angular_rejects_missing_derivative_reference(self):
        with self.assertRaisesRegex(FileNotFoundError, "registration requires"):
            main(["angular", "-C", str(self.path), "--modes", "1", "--angle-step", "90"])
        self.assertFalse((self.path / "angular/angular_intensity_mode1_0.00eV.csv").exists())

    def test_angular_runs_with_registered_reference_without_symmetry(self):
        (self.path / "ref_poscar.vasp").write_text("test\n1\n1 0 0\n0 1 0\n0 0 1\nB\n1\nDirect\n0 0 0\n")
        result = main(["angular", "-C", str(self.path), "--modes", "1", "--angle-step", "90"])
        self.assertEqual(result, 0)
        output = self.path / "angular"
        with (output / "angular_intensity_mode1_0.00eV.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 5)
        self.assertEqual(float(rows[0]["original_parallel"]), 1)
        self.assertEqual(float(rows[0]["original_crossed"]), 0)
        self.assertNotIn("projected_parallel", rows[0])
        self.assertFalse(list(output.glob("symmetry_residuals*")))

    def test_check_exports_only_residuals(self):
        (self.path / "CONTCAR").write_text("test\n1\n1 0 0\n0 1 0\n0 0 1\nB\n1\nDirect\n0 0 0\n")
        symmetry = {"space_group_operations": [{"rotation": np.eye(3).tolist(), "translation": [0, 0, 0]}]}
        (self.path / "symmetry").write_text(yaml.safe_dump(symmetry))
        result = main(["check", "-C", str(self.path)])
        self.assertEqual(result, 0)
        output = self.path / "check"
        with (output / "symmetry_residuals_0.00eV.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(float(row["absolute_residual"]) == 0 for row in rows))
        self.assertFalse(list(output.glob("angular*")))

    def test_complex_eigenvectors_are_opt_in(self):
        data = yaml.safe_load((self.path / "band.yaml").read_text())
        data["phonon"][0]["band"][0]["eigenvector"][0][0] = [0, 1]
        (self.path / "band.yaml").write_text(yaml.safe_dump(data))
        self.assertEqual(read_gamma_modes(self.path / "band.yaml").eigendisplacements[0, 0, 0], 0)
        self.assertEqual(read_gamma_modes(self.path / "band.yaml", complex_eigenvectors=True).eigendisplacements[0, 0, 0], 1j)

    def test_direct_mode_check_is_separate(self):
        (self.path / "CONTCAR").write_text("test\n1\n1 0 0\n0 1 0\n0 0 1\nB\n1\nDirect\n0 0 0\n")
        (self.path / "symmetry").write_text(yaml.safe_dump({"space_group_operations": [
            {"rotation": np.eye(3).tolist(), "translation": [0, 0, 0]}]}))
        (self.path / "irreps.yaml").write_text(yaml.safe_dump({"q-position": [0, 0, 0],
            "normal_modes": [{"band_indices": [i + 1], "ir_label": "A"} for i in range(3)]}))
        self.assertEqual(main(["check", "--kind", "mode", "-C", str(self.path)]), 0)
        output = self.path / "check"
        with (output / "mode_covariance_residuals_0.00eV.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertTrue(all(float(row["absolute_residual"]) == 0 for row in rows))
        self.assertFalse(list(output.glob("symmetry_residuals*")))

    def test_diagnostics_reject_non_gamma_modes(self):
        data = yaml.safe_load((self.path / "band.yaml").read_text())
        data["phonon"][0]["q-position"] = [0.1, 0, 0]
        (self.path / "band.yaml").write_text(yaml.safe_dump(data))
        with self.assertRaises(ValueError):
            read_gamma_modes(self.path / "band.yaml", complex_eigenvectors=True)


if __name__ == "__main__":
    unittest.main()
