"""Helper functions used by nalp_hough.py.

These reproduce the MATLAB utilities referenced (but not shipped) by the
original nalp_hough.m:

    lgwt                  - Legendre-Gauss nodes and weights on [a, b]
    pmn_polynomial_value  - normalized associated Legendre functions
    central_diff          - central-difference derivative on a 1-D grid

The pmn_polynomial_value routine follows John Burkardt's implementation of
the same name.
"""

import numpy as np


def lgwt(N, a, b):
    """Legendre-Gauss nodes x and weights w on the interval [a, b].

    Returns nodes in descending order, matching the MATLAB lgwt convention.
    """
    x, w = np.polynomial.legendre.leggauss(N)  # nodes on [-1, 1], ascending
    # map from [-1, 1] to [a, b]
    x = 0.5 * (a + b) + 0.5 * (b - a) * x
    w = 0.5 * (b - a) * w
    # lgwt returns nodes in descending order
    order = np.argsort(-x)
    return x[order], w[order]


def pmn_polynomial_value(mm, n, m, x):
    """Normalized associated Legendre functions P^m_j(x), j = 0..n.

    Parameters
    ----------
    mm : int
        Number of evaluation points.
    n : int
        Maximum degree.
    m : int
        Order.
    x : (mm,) array_like
        Evaluation points.

    Returns
    -------
    cx : (mm, n + 1) ndarray
        cx[:, j] holds the normalized P^m_j evaluated at x.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    cx = np.zeros((mm, n + 1))

    if m <= n:
        cx[:, m] = 1.0
        factor = 1.0
        for _ in range(1, m + 1):
            cx[:, m] = -factor * cx[:, m] * np.sqrt(1.0 - x ** 2)
            factor += 2.0

    if m + 1 <= n:
        cx[:, m + 1] = x * (2 * m + 1) * cx[:, m]

    for j in range(m + 2, n + 1):
        cx[:, j] = ((2 * j - 1) * x * cx[:, j - 1]
                    + (-j - m + 1) * cx[:, j - 2]) / (j - m)

    # normalization
    for j in range(m, n + 1):
        factor = np.sqrt(((2 * j + 1) * _factorial(j - m))
                         / (2.0 * _factorial(j + m)))
        cx[:, j] = cx[:, j] * factor

    return cx


def _factorial(k):
    """Factorial as a float (values here fit comfortably in double)."""
    result = 1.0
    for i in range(2, k + 1):
        result *= i
    return result


def central_diff(f, x):
    """Central-difference derivative of columns of f with respect to x.

    f may be a 1-D array or a 2-D array whose rows correspond to x.
    Uses second-order central differences in the interior and one-sided
    differences at the endpoints (np.gradient), supporting a non-uniform x.
    """
    f = np.asarray(f, dtype=float)
    x = np.asarray(x, dtype=float).reshape(-1)
    if f.ndim == 1:
        return np.gradient(f, x)
    return np.gradient(f, x, axis=0)
