"""Visual-weight regression checks for absolute Neumann residual plots."""
import unittest
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spectropy_checks import plot_absolute_residual


class ResidualPlotTests(unittest.TestCase):
    def test_intensity_controls_size_and_opacity_and_ticks(self):
        fig, ax = plt.subplots()
        try:
            dots = plot_absolute_residual(ax, [10., 20., 30.], [1., 2., 3.], [0., 1., 4.], 'blue')
            np.testing.assert_allclose(dots.get_sizes(), [5., 27.5, 50.])
            np.testing.assert_allclose(dots.get_facecolors()[:, 3], [1., 1., 1.])
            self.assertEqual(len(ax.collections), 2)  # Dots plus short horizontal ticks, no magnitude series.
            for segment, x, y in zip(ax.collections[1].get_segments(), [10, 20, 30], [1, 2, 3]):
                np.testing.assert_allclose(segment, [[x - 3, y], [x + 3, y]])
        finally:
            plt.close(fig)

    def test_invalid_residual_omitted_and_zero_intensity_visible(self):
        fig, ax = plt.subplots()
        try:
            dots = plot_absolute_residual(ax, [1., 2.], [0., np.nan], [0., 0.])
            self.assertEqual(len(dots.get_offsets()), 1)
            self.assertGreater(dots.get_sizes()[0], 0)
            self.assertGreater(dots.get_facecolors()[0, 3], 0)
            with self.assertRaises(ValueError):
                plot_absolute_residual(ax, [1.], [1.], [-1.])
        finally:
            plt.close(fig)


if __name__ == '__main__':
    unittest.main()
