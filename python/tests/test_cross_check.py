"""Cross-check the two independent Hough solvers against each other.

The Chebyshev and normalized-ALP methods solve the same eigenproblem by
different means, so their physical equivalent depths must agree. This is
the central validation of Wang, Boyd & Akmaev (2016).
"""

import os
import sys

import numpy as np

# Allow running this file directly (python tests/test_cross_check.py): put the
# python/ dir on sys.path so `hough` imports without installing the package.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hough import cheb_hough, nalp_hough  # noqa: E402


def _physical_depths(h, k=7, spurious=1.0e6):
    """Top-k finite, non-spurious equivalent depths, descending.

    The Chebyshev method produces one spurious near-infinite eigenvalue
    (lambda ~ 0), which is filtered out via the ``spurious`` cutoff.
    """
    h = np.asarray(h)
    h = h[np.isfinite(h) & (np.abs(h) < spurious)]
    return np.sort(h)[::-1][:k]


def test_methods_agree_on_equivalent_depths():
    cheb = cheb_hough.compute(s=1.0, sigma=0.5)
    nalp = nalp_hough.compute(s=1, sigma=0.5)

    d_cheb = _physical_depths(cheb.h)
    d_nalp = _physical_depths(nalp.h)

    assert np.allclose(d_cheb, d_nalp, rtol=1e-3), (
        f"\n cheb: {d_cheb}\n nalp: {d_nalp}")


def test_gravest_dw1_depth():
    # The gravest DW1 symmetric mode has an equivalent depth near 0.694 km.
    cheb = cheb_hough.compute(s=1.0, sigma=0.5)
    gravest = _physical_depths(cheb.h, k=1)[0]
    assert abs(gravest - 0.694) < 1e-3, gravest


def _gravest_sw2_wind_peak(res):
    """Peak |hough_u| of the gravest positive-depth (SW2 [2,2]) mode."""
    h = np.where(np.isfinite(res.h) & (np.abs(res.h) < 1.0e6), res.h, np.nan)
    k = int(np.nanargmax(h))
    return np.max(np.abs(res.hough_u[:, k]))


def test_cheb_scalar_is_l2_normalized():
    """cheb_hough L2-normalizes each mode: int hough^2 dx = 1 over [-1,1]
    (Chebyshev-Gauss quadrature), so it shares nalp's physical amplitude."""
    r = cheb_hough.compute(s=2, sigma=1.0)
    n = len(r.x)
    w = (np.pi / n) * np.sqrt(1.0 - r.x ** 2)     # Chebyshev-Gauss weights
    norms = w @ (r.hough ** 2)
    good = np.abs(r.h) < 1.0e6                     # skip spurious modes
    assert np.allclose(norms[good], 1.0, atol=1e-6), norms[good]


def test_cheb_and_nalp_wind_amplitudes_agree():
    """With cheb L2-normalized, the two solvers' wind amplitudes match: the
    gravest SW2 zonal-wind peak is ~2.82 in both (guards the normalization)."""
    peak_cheb = _gravest_sw2_wind_peak(cheb_hough.compute(s=2, sigma=1.0))
    peak_nalp = _gravest_sw2_wind_peak(
        nalp_hough.solve_parity(s=2, sigma=1.0).symmetric)
    assert abs(peak_cheb - peak_nalp) / peak_nalp < 0.02, (peak_cheb, peak_nalp)


if __name__ == "__main__":
    test_methods_agree_on_equivalent_depths()
    test_gravest_dw1_depth()
    test_cheb_scalar_is_l2_normalized()
    test_cheb_and_nalp_wind_amplitudes_agree()
    print("cross-check tests passed")
