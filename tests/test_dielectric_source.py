from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np
from spectropy_config import read_settings, write_settings
from calculate_dielectric_derivatives import read_dielectric, read_dielectrics_xml

def fixture(fields=None, values=None):
    fields=fields or ['energy','xx','yy','zz','xy','yz','zx','yx','zy','xz']
    values=values or [0,1.12345678,2,3,4,5,6,7,8,9]
    array='<array>'+''.join(f'<field>{f}</field>' for f in fields)+'<set><r>'+' '.join(map(str,values))+'</r></set></array>'
    wrong=array.replace('1.12345678','99.0')
    return '<modeling><calculation><dielectricfunction comment="current-current"><real>'+wrong+'</real><imag>'+wrong+'</imag></dielectricfunction><dielectricfunction comment="density-density"><real>'+array+'</real><imag>'+array+'</imag></dielectricfunction></calculation></modeling>'

class DielectricSourceTests(unittest.TestCase):
    def test_settings_xml_only(self):
        with TemporaryDirectory() as d:
            p=Path(d)/'input'
            self.assertEqual(read_settings(p).dielectric_source,'vasprun.xml')
            p.write_text('1 0 0\n1 0 0\nz\ndielectric_source: vasprun.xml\n')
            settings=read_settings(p);write_settings(p,settings)
            self.assertEqual(read_settings(p),settings)
            p.write_text('dielectric_source: outcar\n')
            with self.assertRaisesRegex(ValueError,'must be vasprun.xml'):read_settings(p)

    def test_all_nine_values_are_preserved(self):
        with TemporaryDirectory() as d:
            p=Path(d)/'vasprun.xml';p.write_text(fixture())
            real,imag=read_dielectric(p,0)
            expected=[[1.12345678,4,9],[7,2,5],[6,8,3]]
            np.testing.assert_array_equal(real,expected)
            np.testing.assert_array_equal(imag,expected)
            self.assertNotEqual(real[0,1],real[1,0])
            self.assertNotEqual(real[0,2],real[2,0])
            self.assertNotEqual(real[1,2],real[2,1])
            with self.assertRaisesRegex(ValueError,'readers were removed'):read_dielectric(p,0,'outcar')

    def test_reordered_fields_and_missing_components(self):
        with TemporaryDirectory() as d:
            p=Path(d)/'vasprun.xml'
            p.write_text(fixture(['energy','xz','zy','yx','zx','yz','xy','zz','yy','xx'],[0,9,8,7,6,5,4,3,2,1.12345678]))
            real,_=read_dielectric(p,0)
            np.testing.assert_array_equal(real,[[1.12345678,4,9],[7,2,5],[6,8,3]])
            p.write_text(fixture(['energy','xx','yy','zz','xy','yz','zx'],[0,1,2,3,4,5,6]))
            with self.assertRaisesRegex(ValueError,'all nine'):read_dielectric(p,0)

    def test_invalid_and_out_of_grid(self):
        with TemporaryDirectory() as d:
            p=Path(d)/'vasprun.xml';p.write_text(fixture())
            with self.assertRaisesRegex(ValueError,'outside'):read_dielectrics_xml(p,[1])
            p.write_text('<modeling/>')
            with self.assertRaisesRegex(ValueError,'density-density'):read_dielectric(p,0)

if __name__=='__main__':unittest.main()
