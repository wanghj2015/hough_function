"""Cross-check the two independent Hough solvers against each other.

The Chebyshev and normalized-ALP methods solve the same eigenproblem by
different means, so their physical equivalent depths must agree. This is
the central validation of Wang, Boyd & Akmaev (2016).
"""

import numpy as np

from hough import cheb_hough, nalp_hough


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


if __name__ == "__main__":
    test_methods_agree_on_equivalent_depths()
    test_gravest_dw1_depth()
    print("cross-check tests passed")
