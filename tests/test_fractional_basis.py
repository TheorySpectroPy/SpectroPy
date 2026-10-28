import copy
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

import reconstruct_dielectric_derivatives as reconstruction
from spectropy_displacements import Displacement
import spectropy_displacements as generation
import spectropy_derivatives as derivatives
from spectropy_structure import Structure


class SavedFractionalBasisTests(unittest.TestCase):
    def setUp(self):
        self.lattice = np.array([[1., 0., 0.], [-0.5, np.sqrt(3) / 2, 0.], [0., 0., 2.]])
        rotation = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
        self.rotations = np.array([np.eye(3), rotation, rotation @ rotation])
        self.symmetry = {1: {"rotations": self.rotations}}

    def check_equivalence(self, vectors, symmetry):
        displacements = [Displacement(1, np.asarray(v) * 0.01, chr(97 + i)) for i, v in enumerate(vectors)]
        bases = reconstruction.build_fractional_bases(displacements, self.lattice, symmetry)
        rng = np.random.default_rng(42)
        values = [(d.fractional_vector @ self.lattice, rng.normal(size=(3, 3)), rng.normal(size=(3, 3))) for d in displacements]
        args = ({1: values}, self.lattice, symmetry, {1: 1}, {})
        legacy = reconstruction.reconstruct_fractional_derivatives(*args)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fractional_basis.json"
            reconstruction.write_fractional_bases(path, bases)
            saved = reconstruction.read_fractional_bases(path, displacements, self.lattice, symmetry)
            with patch.object(reconstruction, "expand_direction_orbit", side_effect=AssertionError("Repeated orbit search")), \
                 patch.object(reconstruction, "_basis_indices", side_effect=AssertionError("Repeated basis selection")):
                result = reconstruction.reconstruct_fractional_derivatives(*args, bases=saved)
            for old, new in zip(legacy, result):
                np.testing.assert_array_equal(old[1], new[1])
            stale_lattice = self.lattice.copy()
            stale_lattice[2, 2] += 0.01
            with self.assertRaisesRegex(ValueError, "stale"):
                reconstruction.read_fractional_bases(path, displacements, stale_lattice, symmetry)
            changed = [Displacement(1, displacements[0].fractional_vector * 2, "a"), *displacements[1:]]
            with self.assertRaisesRegex(ValueError, "stale"):
                reconstruction.read_fractional_bases(path, changed, self.lattice, symmetry)
            changed_symmetry = copy.deepcopy(symmetry)
            changed_symmetry[1]["rotations"] = np.eye(3)[None]
            if not np.array_equal(symmetry[1]["rotations"], changed_symmetry[1]["rotations"]):
                with self.assertRaisesRegex(ValueError, "stale"):
                    reconstruction.read_fractional_bases(path, displacements, self.lattice, changed_symmetry)

    def test_one_mixed_seed(self):
        self.check_equivalence([[1, 0, 1], [-1, 0, -1]], self.symmetry)

    def test_two_seeds_and_synthesized_sign(self):
        inversion_z = np.diag([1, 1, -1])
        rotations = np.concatenate([self.rotations, [r @ inversion_z for r in self.rotations]])
        self.check_equivalence([[1, 0, 0], [-1, 0, 0], [0, 0, 1]], {1: {"rotations": rotations}})

    def test_three_seeds_trivial_site(self):
        self.check_equivalence([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], {1: {"rotations": np.eye(3)[None]}})

    def test_incomplete_basis_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "span three dimensions"):
            reconstruction.build_fractional_bases([Displacement(1, np.array([0.01, 0, 0]), "a")], self.lattice, {1: {"rotations": np.eye(3)[None]}})

    def test_generation_writes_reusable_basis(self):
        structure = Structure("test", self.lattice, ("B",), np.array([1]), np.zeros((1, 3)))
        previous = Path.cwd()
        with tempfile.TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                with patch.object(generation, "read_structure", return_value=structure), \
                     patch.object(generation, "_amplitude", return_value=np.array([0.03] * 3)), \
                     patch.object(generation, "read_symmetry_file", return_value=(None, None, {1: 1})), \
                     patch.object(generation, "read_site_symmetries", return_value=self.symmetry):
                    generation.run_displacements("fractional")
                records = generation.read_displacements()
                bases = reconstruction.read_fractional_bases("fractional_basis.json", records, self.lattice, self.symmetry)
                self.assertEqual(len(bases["atoms"]["1"]["entries"]), 3)
                self.assertEqual(len(records), 2)
                Path("CONTCAR").touch()
                Path("symmetry").touch()
                measured = {1: [(d.fractional_vector @ self.lattice, np.eye(3), np.eye(3) * 2) for d in records]}
                with patch.object(derivatives, "read_structure", return_value=structure), \
                     patch.object(derivatives, "read_symmetry_file", return_value=(self.rotations, np.zeros((3, 3)), {1: 1})), \
                     patch.object(derivatives, "read_site_symmetries", return_value=self.symmetry), \
                     patch.object(derivatives, "compute_mapping_matrices", return_value=(None, None, {})), \
                     patch.object(derivatives, "collect_vasprun_files"), \
                     patch.object(derivatives, "read_settings", return_value=SimpleNamespace(laser_energies=[0., 1.5], difference_scheme='central', dielectric_source='outcar')), \
                     patch.object(derivatives, "read_displaced_dielectrics", return_value=measured), \
                     patch.object(reconstruction, "expand_direction_orbit", side_effect=AssertionError("Repeated orbit search")):
                    derivatives.run_derivatives(use_saved_basis=True)
                self.assertTrue(Path("epsilon_derivative_0.00").is_file())
                self.assertTrue(Path("epsilon_derivative_1.50").is_file())
                # Default manual path must work even with an unusable lookup.
                Path("fractional_basis.json").write_text("not valid json")
                with patch.object(derivatives, "read_structure", return_value=structure), \
                     patch.object(derivatives, "read_symmetry_file", return_value=(self.rotations, np.zeros((3, 3)), {1: 1})), \
                     patch.object(derivatives, "read_site_symmetries", return_value=self.symmetry), \
                     patch.object(derivatives, "compute_mapping_matrices", return_value=(None, None, {})), \
                     patch.object(derivatives, "collect_vasprun_files"), \
                     patch.object(derivatives, "read_settings", return_value=SimpleNamespace(laser_energies=[0.], difference_scheme='central', dielectric_source='outcar')), \
                     patch.object(derivatives, "read_displaced_dielectrics", return_value=measured):
                    derivatives.run_derivatives()
            finally:
                os.chdir(previous)


if __name__ == "__main__":
    unittest.main()
