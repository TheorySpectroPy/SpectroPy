"""Raw spectra must survive optional site projection and multiple energies."""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import calculate_spectrum as spectrum
from spectropy_config import Settings, write_settings


class RawSpectrumOutputsTests(unittest.TestCase):
    def run_case(self, denoise, legacy=False):
        settings = Settings(incident_polarization=(1, 0, 0),
                            scattered_polarization=(1, 0, 0), surface_normal='z',
                            site_symmetry_projection=denoise, laser_energies=(0., 1.96),
                            polarization='average')
        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            try:
                os.chdir(tmp)
                write_settings('input', settings)
                for energy in settings.laser_energies:
                    name = f'epsilon_derivative_{energy:.2f}'
                    Path(name if legacy else 'raw_' + name).touch()
                    # Stale corrected files must not activate projection.
                    Path('corrected_' + name).touch()
                tensors = np.zeros((1, 3, 3, 3))
                tensors[0, 0, 0, 0] = 2.

                def derivatives(path, natoms):
                    scale = 0.5 if str(path).startswith('corrected_') else 1.
                    return scale * tensors, np.zeros_like(tensors)

                with patch.object(spectrum, 'read_band_yaml', return_value=(
                        np.ones(3), np.eye(3).reshape(3, 1, 3), 1)), \
                     patch.object(spectrum, 'read_irreps', return_value=None), \
                     patch.object(spectrum, 'read_dielectric_derivatives', side_effect=derivatives):
                    spectrum.run_raman_tensors_for_input()
                for energy in settings.laser_energies:
                    suffix = f'{energy:.2f}eV'
                    raw = Path(f'raw_Raman_tensor_{suffix}').read_text()
                    self.assertTrue(Path(f'raw_Raman_intensity_polarization_averaged_{suffix}').exists())
                    self.assertIn('2.000', raw)
                    active = Path(f'Raman_tensor_{suffix}').read_text()
                    if denoise:
                        corrected = Path(f'corrected_Raman_tensor_{suffix}').read_text()
                        self.assertNotEqual(raw, corrected)
                        self.assertEqual(active, corrected)
                        raw_i = np.loadtxt(f'raw_Raman_intensity_complex_{suffix}', usecols=1)
                        clean_i = np.loadtxt(f'corrected_Raman_intensity_complex_{suffix}', usecols=1)
                        np.testing.assert_allclose(raw_i, 4 * clean_i, atol=5e-4)
                    else:
                        self.assertEqual(active, raw)
                        self.assertFalse(Path(f'corrected_Raman_tensor_{suffix}').exists())
            finally:
                os.chdir(previous)

    def test_enabled_preserves_raw_and_corrected(self):
        self.run_case(True)

    def test_disabled_produces_raw_only(self):
        self.run_case(False)

    def test_legacy_unprojected_input(self):
        self.run_case(False, legacy=True)

    def test_enabled_rejects_missing_raw(self):
        with self.assertRaisesRegex(FileNotFoundError, 'raw_epsilon_derivative'):
            self.run_case(True, legacy=True)


if __name__ == '__main__':
    unittest.main()
