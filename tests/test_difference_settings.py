"""Input-file selection of finite differences, with CLI/API overrides."""
from pathlib import Path
import os
import numpy as np
import tempfile
import unittest
from unittest.mock import patch
from spectropy_config import read_settings, write_settings
from spectropy_derivatives import run_derivatives
from spectropy_displacements import run_displacements, read_displacements
from spectropy_structure import read_structure


class DifferenceSettingsTests(unittest.TestCase):
    def test_default_and_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input'
            self.assertEqual(read_settings(path).difference_scheme, 'central')
            path.write_text('1 0 0\n1 0 0\nz\ndifference_scheme: forward\n')
            settings = read_settings(path)
            self.assertEqual(settings.difference_scheme, 'forward')
            self.assertFalse(hasattr(settings, 'equilibrium_vasprun'))
            write_settings(path, settings)
            self.assertEqual(read_settings(path), settings)

    def test_invalid_scheme(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input'
            path.write_text('difference_scheme: backwards\n')
            with self.assertRaisesRegex(ValueError, 'central or forward'):
                read_settings(path)

    def test_derivatives_obey_input_and_explicit_override(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input'
            path.write_text('difference_scheme: forward\n')
            with self.assertRaisesRegex(ValueError, 'ra_pos_equilibrium'):
                run_derivatives(path)
            # Central override must pass the forward-only reference validation.
            with patch('spectropy_derivatives.os.path.exists', return_value=False):
                with self.assertRaises(SystemExit):
                    run_derivatives(path, difference_scheme='central')

    def test_generated_reference_and_automatic_derivative_lookup(self):
        original = Path.cwd()
        for scheme in ('central', 'forward'):
            with tempfile.TemporaryDirectory() as directory:
                try:
                    os.chdir(directory)
                    Path('CONTCAR').write_text('test\n1\n2 0 0\n0 2 0\n0 0 2\nB\n1\nDirect\n0.2 0.3 0.4\n')
                    Path('input').write_text(f'1 0 0\n1 0 0\nz\ndisplacement_amplitude: 0.01\ndifference_scheme: {scheme}\n')
                    records = run_displacements(mode='full')
                    self.assertEqual(len(read_displacements()), len(records))
                    self.assertTrue(all(np.linalg.norm(r.fractional_vector) > 0 for r in records))
                    self.assertEqual(Path('displacements/ra_pos_equilibrium').exists(), scheme == 'forward')
                    if scheme == 'forward':
                        np.testing.assert_allclose(read_structure('displacements/ra_pos_equilibrium/POSCAR').fractional_positions,
                                                   read_structure('CONTCAR').fractional_positions)
                        Path('displacements/ra_pos_equilibrium/OUTCAR').touch()
                        measurements = {1: [(r.fractional_vector * 2, np.eye(3), np.zeros((3, 3))) for r in records]}
                        with patch('spectropy_derivatives.collect_vasprun_files'), \
                             patch('spectropy_derivatives.read_dielectric', return_value=(np.eye(3), np.zeros((3, 3)))) as reader, \
                             patch('spectropy_derivatives.read_displaced_dielectrics', return_value=measurements):
                            run_derivatives()
                        reader.assert_called_once_with('displacements/ra_pos_equilibrium/OUTCAR', 0., 'outcar')
                        self.assertTrue(Path('epsilon_derivative_0.00').exists())
                finally:
                    os.chdir(original)


if __name__ == '__main__':
    unittest.main()
