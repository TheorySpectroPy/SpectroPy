import unittest
import numpy as np
from reconstruct_dielectric_derivatives import compute_atom_tensor, compute_fractional_atom_tensor


class ForwardDifferenceTests(unittest.TestCase):
    def setUp(self):
        self.eq = np.eye(3) * 2
        self.derivative = np.arange(27.).reshape(3, 3, 3)
        self.identity = np.eye(3)[None]

    def test_cartesian_linear_and_negative_only(self):
        for sign in (1, -1):
            measured = [(sign * .01 * np.eye(3)[i], self.eq + sign * .01 * self.derivative[:, :, i]) for i in range(3)]
            actual = compute_atom_tensor(measured, self.identity, difference_scheme="forward", equilibrium=self.eq)
            np.testing.assert_allclose(actual, self.derivative, atol=1e-12)

    def test_fractional_linear_skew_lattice(self):
        lattice = np.array([[2., 0., 0.], [-1., 1.7, 0.], [0., 0., 5.]])
        measured = []
        for i in range(3):
            displacement = .01 * np.eye(3)[i] @ lattice
            measured.append((displacement, self.eq + np.einsum('abi,i->ab', self.derivative, displacement)))
        actual = compute_fractional_atom_tensor(measured, self.identity, lattice, difference_scheme="forward", equilibrium=self.eq)
        np.testing.assert_allclose(actual, self.derivative, atol=1e-12)

    def test_forward_first_order_central_second_order(self):
        errors = []
        for step in (.02, .01):
            measured = [(sign * step * np.eye(3)[i], self.eq + sign * step * self.derivative[:, :, i] + step**2) for i in range(3) for sign in (1, -1)]
            central = compute_atom_tensor(measured, self.identity)
            forward = compute_atom_tensor(measured, self.identity, difference_scheme="forward", equilibrium=self.eq)
            np.testing.assert_allclose(central, self.derivative, atol=1e-12)
            errors.append(np.linalg.norm(forward - self.derivative))
        self.assertAlmostEqual(errors[0] / errors[1], 2)

    def test_missing_equilibrium_rejected(self):
        with self.assertRaisesRegex(ValueError, "equilibrium"):
            compute_atom_tensor([], self.identity, difference_scheme="forward")


if __name__ == '__main__':
    unittest.main()
