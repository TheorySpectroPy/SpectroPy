import numpy as np
import pytest
from generate_raman_plots import broaden_spectrum


def test_narrow_line_grid_resolves_wide_spectral_range():
    grid, values = broaden_spectrum(np.array([.0, 799.713, 1600.]), np.array([0.,1.,0.]), 1., 'l')
    assert np.max(np.diff(grid)) <= .05+1e-10
    assert values.max() >= 1/(1+.05**2)-1e-10


def test_invalid_linewidth_rejected():
    for width in [0., -1., float('nan')]:
        with pytest.raises(ValueError):
            broaden_spectrum(np.array([1.]),np.array([1.]),width,'l')
