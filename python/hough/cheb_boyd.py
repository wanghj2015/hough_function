"""cheb_boyd - Compute differential matrices for the Chebyshev collocation
method, with an optional parity factor (pf).

Python port of cheb_boyd.m from:
Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. Geoscientific Model Development, 9(4), 1477-1488.
"""

import numpy as np


def cheb_boyd(N, pf):
    """Return first (D1) and second (D2) derivative matrices and the grid x.

    Parameters
    ----------
    N : int
        Number of collocation points.
    pf : int
        Parity factor (0 or 1).

    Returns
    -------
    D1, D2 : (N, N) ndarray
        First- and second-derivative collocation matrices.
    x : (N,) ndarray
        Collocation points.
    """
    t = (np.pi / (2 * N) * np.arange(1, 2 * N, 2)).reshape(-1, 1)  # (N, 1)
    x = np.cos(t).ravel()
    n = np.arange(N)

    ss = np.sin(t)
    cc = np.cos(t)
    sx = np.tile(ss, (1, N))          # (N, N)
    cx = np.tile(cc, (1, N))
    nx = np.tile(n, (N, 1))           # (N, N)
    tx = np.tile(t, (1, N))
    tn = np.cos(nx * tx)

    if pf == 0:
        phi2 = tn
        PT = -nx * np.sin(nx * tx)
        phiD2 = -PT / sx
        PTT = -nx ** 2 * tn
        phiDD2 = (sx * PTT - cx * PT) / sx ** 3
    else:
        phi2 = tn * sx
        PT = -nx * np.sin(nx * tx) * sx + tn * cx
        phiD2 = -PT / sx
        PTT = (-nx ** 2 * tn * sx
               - 2 * nx * np.sin(nx * tx) * cx - tn * sx)
        phiDD2 = (sx * PTT - cx * PT) / sx ** 3

    # MATLAB "A / phi2" is A * inv(phi2), i.e. solve X * phi2 = A.
    D1 = np.linalg.solve(phi2.T, phiD2.T).T   # first derivatives
    D2 = np.linalg.solve(phi2.T, phiDD2.T).T  # second derivatives
    return D1, D2, x
