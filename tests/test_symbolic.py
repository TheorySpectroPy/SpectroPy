import unittest

import sympy as sp

from spectropy_symbolic import derive_constraints, exact_rotation, tensor_action


class SymbolicConstraintsTests(unittest.TestCase):
    def test_identity_and_optional_reciprocity(self):
        ops = [([0], sp.eye(3))]
        self.assertEqual(derive_constraints(1, ops)['parameter_count'], 27)
        self.assertEqual(derive_constraints(1, ops, True)['parameter_count'], 18)

    def test_inversion_center_forbids_polar_rank_three(self):
        result = derive_constraints(1, [([0], sp.eye(3)), ([0], -sp.eye(3))])
        self.assertEqual(result['parameter_count'], 0)
        self.assertEqual(result['forms'][0], sp.zeros(27, 1))

    def test_inversion_exchanges_atoms_not_zero(self):
        result = derive_constraints(2, [([0, 1], sp.eye(3)), ([1, 0], -sp.eye(3))], True)
        self.assertEqual(result['parameter_count'], 18)
        self.assertEqual(result['forms'][1], -result['forms'][0])

    def test_horizontal_mirror_three_index_parity(self):
        result = derive_constraints(1, [([0], sp.eye(3)), ([0], sp.diag(1, 1, -1))], True)
        self.assertEqual(result['parameter_count'], 10)
        v = result['forms'][0]
        for i in range(3):
            for a in range(3):
                for b in range(3):
                    if sum(k == 2 for k in (i, a, b)) % 2:
                        self.assertEqual(v[9*i+3*a+b], 0)
        # Out-of-plane displacement CAN have xz/yz response despite this mirror.
        self.assertNotEqual(v[18+2], 0)

    def test_action_index_order_matches_explicit_contraction(self):
        s = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        v = sp.Matrix(range(27))
        actual = tensor_action(s)*v
        for i in range(3):
            for a in range(3):
                for b in range(3):
                    expected = sum(s[i,j]*s[a,c]*s[b,d]*v[9*j+3*c+d]
                                   for j in range(3) for c in range(3) for d in range(3))
                    self.assertEqual(actual[9*i+3*a+b], expected)

    def test_rejects_bad_rotation_or_incomplete_group(self):
        with self.assertRaises(ValueError):
            exact_rotation(sp.diag(1, 1, 2))
        c4 = sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
        with self.assertRaisesRegex(ValueError, 'closed'):
            derive_constraints(1, [([0], sp.eye(3)), ([0], c4)])

    def test_mos2_site_counts_and_sulfur_mirror_transfer(self):
        c3 = sp.Matrix([[-sp.Rational(1,2), -sp.sqrt(3)/2, 0],
                        [sp.sqrt(3)/2, -sp.Rational(1,2), 0], [0,0,1]])
        vertical = sp.diag(-1, 1, 1)
        horizontal = sp.diag(1, 1, -1)
        ops = []
        for power in range(3):
            for mirror in range(2):
                s = c3**power * vertical**mirror
                ops.append(([0,1,2], s))
                ops.append(([0,2,1], horizontal*s))
        result = derive_constraints(3, ops, True)
        self.assertEqual([len(o['parameters']) for o in result['orbits']], [1,4])
        self.assertEqual(result['parameter_count'], 5)
        self.assertEqual(result['forms'][0][18:27], [0]*9)
        self.assertEqual(result['forms'][2], tensor_action(horizontal)*result['forms'][1])
        self.assertNotEqual(result['forms'][1][2], 0)  # legitimate S dx -> epsilon_xz

    def test_floating_crystallographic_rotation_recovery(self):
        exact = sp.Matrix([[-sp.Rational(1,2), -sp.sqrt(3)/2, 0],
                           [sp.sqrt(3)/2, -sp.Rational(1,2), 0], [0,0,1]])
        floats = [[float(x) for x in exact.row(i)] for i in range(3)]
        self.assertEqual(exact_rotation(floats), exact)


if __name__ == '__main__':
    unittest.main()
