from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import os
import unittest
import xml.etree.ElementTree as ET
import numpy as np

import calculate_dielectric_derivatives as reader
import spectropy_derivatives as workflow
from spectropy_displacements import Displacement


def xml_fixture():
    fields = ['energy', 'xz', 'zy', 'yx', 'zx', 'yz', 'xy', 'zz', 'yy', 'xx']
    def array(offset):
        rows = ''.join('<r>'+' '.join(map(str, [e]+[offset+10*e+i for i in range(1,10)]))+'</r>'
                       for e in [0., 1., 2.])
        return '<array>'+''.join(f'<field>{f}</field>' for f in fields)+'<set>'+rows+'</set></array>'
    return ('<modeling><dielectricfunction comment="density-density"><real>'+array(0)+
            '</real><imag>'+array(100)+'</imag></dielectricfunction></modeling>')


class BatchedDielectricTests(unittest.TestCase):
    def test_identical_samples_and_one_parse_per_file(self):
        with TemporaryDirectory() as directory:
            p = Path(directory)
            ds = [Displacement(1, np.array([.01,0,0]), 'a'),
                  Displacement(2, np.array([0,.01,0]), 'a')]
            for d in ds:
                (p/f'{d.name}.xml').write_text(xml_fixture())
            energies = [2., .49, 1.5, 0., 2.]
            expected = [reader.read_displaced_dielectrics(ds, np.eye(3), e, p) for e in energies]
            with patch.object(reader.ET, 'parse', wraps=ET.parse) as parse:
                actual = reader.read_displaced_dielectrics_many(ds, np.eye(3), energies, p)
            self.assertEqual(parse.call_count, len(ds))
            for old, new in zip(expected, actual):
                self.assertEqual(list(old), list(new))
                for atom in old:
                    for a, b in zip(old[atom][0], new[atom][0]):
                        np.testing.assert_array_equal(a, b)
                    self.assertFalse(np.allclose(new[atom][0][1], new[atom][0][1].T))
                    self.assertTrue(new[atom][0][1].flags.owndata)
            self.assertFalse(np.shares_memory(actual[0][1][0][1], actual[-1][1][0][1]))
            self.assertEqual(reader.read_displaced_dielectrics_many(ds, np.eye(3), [], p), [])

    def test_validation_is_not_bypassed(self):
        with TemporaryDirectory() as directory:
            p = Path(directory)
            d = Displacement(1, np.array([.01,0,0]), 'a')
            f = p/'atom1a.xml'
            f.write_text(xml_fixture())
            with self.assertRaisesRegex(ValueError, 'outside'):
                reader.read_displaced_dielectrics_many([d], np.eye(3), [0.,3.], p)
            f.write_text(xml_fixture()[:-6])
            with self.assertRaises(ET.ParseError):
                reader.read_displaced_dielectrics_many([d], np.eye(3), [0.], p)
            f.write_text(xml_fixture().replace('<field>zy</field>', '<field>xx</field>'))
            with self.assertRaisesRegex(ValueError, 'all nine'):
                reader.read_displaced_dielectrics_many([d], np.eye(3), [0.], p)
            f.write_text(xml_fixture().replace('101.0', 'nan'))
            with self.assertRaisesRegex(ValueError, 'Invalid'):
                reader.read_displaced_dielectrics_many([d], np.eye(3), [0.], p)
            with self.assertRaisesRegex(ValueError, 'removed'):
                reader.read_displaced_dielectrics_many([d], np.eye(3), [0.], p, source='OUTCAR')

    def test_forward_reference_and_displacements_batched_once(self):
        with TemporaryDirectory() as directory:
            old_cwd = os.getcwd()
            os.chdir(directory)
            try:
                for name in ['CONTCAR','displacements.dat','equilibrium.xml']:
                    Path(name).touch()
                structure = SimpleNamespace(lattice=np.eye(3), fractional_positions=np.zeros((1,3)))
                settings = SimpleNamespace(site_symmetry_projection=False, dielectric_source='vasprun.xml',
                                           difference_scheme='forward', laser_energies=[0.,1.])
                ds = [object()]
                eq = [(np.eye(3),np.zeros((3,3))), (2*np.eye(3),np.zeros((3,3)))]
                measurements = [{1: []}, {1: []}]
                with patch.object(workflow,'read_settings',return_value=settings), \
                     patch.object(workflow,'read_structure',return_value=structure), \
                     patch.object(workflow,'read_displacements',return_value=ds), \
                     patch.object(workflow,'read_displacement_mode',return_value='full'), \
                     patch.object(workflow,'collect_vasprun_files') as collect, \
                     patch.object(workflow,'read_dielectrics_xml',return_value=eq) as reference, \
                     patch.object(workflow,'read_displaced_dielectrics_many',return_value=measurements) as batch, \
                     patch.object(workflow,'_write_dielectric_log'), \
                     patch.object(workflow,'write_derivative_outputs') as write, \
                     patch.object(workflow,'reconstruct_fractional_derivatives',return_value=({},{})) as reconstruct:
                    workflow.run_derivatives(equilibrium_vasprun='equilibrium.xml')
                    collect.assert_called_once()
                    reference.assert_called_once_with('equilibrium.xml', [0.,1.])
                    batch.assert_called_once()
                    self.assertEqual(write.call_count, 2)
                    for call, pair in zip(reconstruct.call_args_list, eq):
                        for a,b in zip(call.kwargs['equilibrium'],pair):
                            np.testing.assert_array_equal(a,b)
            finally:
                os.chdir(old_cwd)


if __name__ == '__main__':
    unittest.main()
