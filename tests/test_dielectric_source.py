from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import numpy as np

from spectropy_config import read_settings, write_settings
from calculate_dielectric_derivatives import read_dielectric, read_diel_from_outcar
from spectropy_derivatives import _write_derivatives
from calculate_spectrum import read_dielectric_derivatives


class DielectricSourceTests(unittest.TestCase):
    def test_settings_default_xml_and_invalid(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'input'
            self.assertEqual(read_settings(path).dielectric_source, 'outcar')
            path.write_text('1 0 0\n1 0 0\nz\ndielectric_source: vasprun.xml\n')
            settings = read_settings(path)
            write_settings(path, settings)
            self.assertEqual(read_settings(path), settings)
            path.write_text('dielectric_source: unknown\n')
            with self.assertRaisesRegex(ValueError, 'outcar or vasprun.xml'):
                read_settings(path)

    def test_precision_components_and_grid(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'OUTCAR'
            path.write_text(' frequency dependent REAL DIELECTRIC FUNCTION density-density\n'
                            'header\nseparator\n0.000000 1.745363 1.745393 1.274160 0.000001 -0.000003 -0.000002\n\n'
                            ' frequency dependent IMAGINARY DIELECTRIC FUNCTION density-density\n'
                            'header\nseparator\n0.000000 0 0 0 0 0 0\n\n'
                            ' frequency dependent REAL DIELECTRIC FUNCTION current-current\n'
                            'header\nseparator\n0.000000 9 9 9 9 9 9\n')
            real, imag = read_dielectric(path, 0., 'outcar')
            np.testing.assert_array_equal(real, [[1.745363, .000001, -.000002],
                                                [.000001, 1.745393, -.000003],
                                                [-.000002, -.000003, 1.274160]])
            self.assertFalse(imag.any())
            with self.assertRaisesRegex(ValueError, 'not present'):
                read_diel_from_outcar(path, .001)
            with patch('calculate_dielectric_derivatives.read_diel_from_xml') as xml:
                with self.assertRaises(FileNotFoundError):
                    read_dielectric(path.parent / 'missing', 0., 'outcar')
                xml.assert_not_called()
                read_dielectric('vasprun.xml', 1.505, 'vasprun.xml')
                xml.assert_called_once_with('vasprun.xml', 1.505)

    def test_derivative_export_retains_precision(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'derivative'
            real = {1: np.arange(27).reshape(3, 3, 3) * .00000123456789012345}
            imag = {1: real[1] / 3}
            _write_derivatives(path, 0., np.eye(3), np.zeros((1, 3)), real, imag)
            a, b = read_dielectric_derivatives(path, 1)
            np.testing.assert_allclose(a[0], real[1].transpose(2, 0, 1)/(4*np.pi), rtol=1e-14, atol=0)
            np.testing.assert_allclose(b[0], imag[1].transpose(2, 0, 1)/(4*np.pi), rtol=1e-14, atol=0)
