from pathlib import Path
import tempfile
import unittest
import numpy as np
from calculate_spectrum import read_dielectric_derivatives


class DerivativeIndexingTests(unittest.TestCase):
    def test_explicit_atom_labels_not_block_order(self):
        for label in ('Atom {atom}', 'atom{atom}'):
            blocks = []
            for part, multiplier in (('Real', 1), ('Imaginary', 10)):
                blocks.append(f'! {part} Part\n')
                for atom in (2, 1):
                    blocks.append(label.format(atom=atom) + '\n')
                    blocks.extend((' '.join([str(atom * multiplier)] * 9) + '\n') for _ in range(3))
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'epsilon_derivative'
                path.write_text(''.join(blocks))
                real, imag = read_dielectric_derivatives(path, 2)
            np.testing.assert_array_equal(real[0], 1)
            np.testing.assert_array_equal(real[1], 2)
            np.testing.assert_array_equal(imag[0], 10)
            np.testing.assert_array_equal(imag[1], 20)


if __name__ == '__main__':
    unittest.main()
