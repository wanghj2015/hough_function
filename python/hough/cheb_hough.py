"""Compute Hough functions using Chebyshev collocation methods.

Python port of cheb_hough.m from:
Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. Geoscientific Model Development, 9(4), 1477-1488.
"""

from collections import namedtuple

import numpy as np

from .cheb_boyd import cheb_boyd

# Physical constants
A_EARTH = 6.370e6                       # Earth radius (m)
G = 9.81                                # gravity (m s^-2)
OMEGA = 2.0 * np.pi / (24.0 * 3600.0)   # Earth rotation rate (s^-1)

HoughResult = namedtuple(
    "HoughResult", ["x", "lamb", "hough", "hough_u", "hough_v", "h"])


def compute(s=1.0, sigma=0.5, N=62):
    """Compute Hough functions via Chebyshev collocation.

    Parameters
    ----------
    s : float
        Zonal wavenumber. Common cases: DW1 (s=1, sigma=0.5),
        SW2 (s=2, sigma=1.0), TW3 (s=3, sigma=1.5).
    sigma : float
        Normalized frequency (frequency / (2*OMEGA)).
    N : int
        Number of collocation points.

    Returns
    -------
    HoughResult
        Fields: x (grid), lamb (eigenvalues), hough (mode shapes,
        columns are modes), hough_u / hough_v (wind components),
        h (equivalent depth, km).
    """
    parity_factor = int(s) % 2
    D1, D2, x = cheb_boyd(N, parity_factor)

    a2 = (1 - x ** 2) / (sigma ** 2 - x ** 2)
    a1 = 2.0 * x * (1 - sigma ** 2) / (sigma ** 2 - x ** 2) ** 2
    a0 = -1.0 / (sigma ** 2 - x ** 2) * (
        (s / sigma) * (sigma ** 2 + x ** 2) / (sigma ** 2 - x ** 2)
        + s ** 2 / (1 - x ** 2))

    A = np.diag(a2) @ D2 + np.diag(a1) @ D1 + np.diag(a0)

    d, v = np.linalg.eig(A)
    lamb = np.real(d)

    # sort eigenvalues and -vectors (descending lambda)
    ii = np.argsort(-lamb)
    lamb = lamb[ii]
    hough = np.real(v[:, ii])

    # equivalent depth (km)
    h = -4.0 * A_EARTH ** 2 * OMEGA ** 2 / G / lamb / 1000.0

    # Hough functions for wind components
    b1 = (sigma ** 2 - x ** 2) * np.sqrt(1.0 - x ** 2)
    b2 = np.sqrt(1.0 - x ** 2) / (sigma ** 2 - x ** 2)
    hough_u = (np.diag(s / b1) @ hough
               - np.diag(b2 * x / sigma) @ D1 @ hough)
    hough_v = (np.diag((s / sigma) * x / b1) @ hough
               - np.diag(b2) @ D1 @ hough)

    return HoughResult(x, lamb, hough, hough_u, hough_v, h)


if __name__ == "__main__":
    res = compute()
    print("first equivalent depths (km):", np.round(res.h[:6], 4))
