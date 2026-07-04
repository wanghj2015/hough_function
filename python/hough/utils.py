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


def jacobi_eigenvalue(a_in, it_max=100):
    """Classical cyclic Jacobi rotation eigensolver for a real symmetric
    matrix, faithfully ported from John Burkardt's jacobi_eigenvalue.m
    (https://people.sc.fsu.edu/~jburkardt/m_src/jacobi_eigenvalue/), the
    algorithm the paper cites for solving the F1/F2 matrices.

    Unlike LAPACK's eig/eigh (used by MATLAB's `eig` and NumPy's `eigh`),
    Jacobi rotations start from V = identity and build up the eigenvectors
    through an explicit, deterministic sequence of plane rotations. That
    gauge is what reproduces the published PDF's curve signs -- verified
    against every mode in Figs. 1-3 by pixel-level digitization of the paper
    -- whereas eig/eigh's gauge only sometimes agrees with it. See
    docs/reference.md.

    Returns (eigenvalues, eigenvectors) with eigenvectors as columns,
    unsorted (same convention as numpy.linalg.eigh).
    """
    n = a_in.shape[0]
    a = a_in.copy()
    v = np.eye(n)
    d = np.diag(a).copy()
    bw = d.copy()
    zw = np.zeros(n)
    it_num = 0

    while it_num < it_max:
        it_num += 1
        thresh = np.sqrt(np.sum(np.triu(a, 1) ** 2)) / (4 * n)
        if thresh == 0.0:
            break

        for p in range(n):
            for q in range(p + 1, n):
                gapq = 10.0 * abs(a[p, q])
                termp = gapq + abs(d[p])
                termq = gapq + abs(d[q])

                if 4 < it_num and termp == abs(d[p]) and termq == abs(d[q]):
                    a[p, q] = 0.0
                elif thresh <= abs(a[p, q]):
                    h = d[q] - d[p]
                    term = abs(h) + gapq
                    if term == abs(h):
                        t = a[p, q] / h
                    else:
                        theta = 0.5 * h / a[p, q]
                        t = 1.0 / (abs(theta) + np.sqrt(1.0 + theta * theta))
                        if theta < 0.0:
                            t = -t
                    c = 1.0 / np.sqrt(1.0 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[p, q]
                    zw[p] -= h
                    zw[q] += h
                    d[p] -= h
                    d[q] += h
                    a[p, q] = 0.0

                    for i in range(p):
                        g, h = a[i, p], a[i, q]
                        a[i, p] = g - s * (h + g * tau)
                        a[i, q] = h + s * (g - h * tau)
                    for i in range(p + 1, q):
                        g, h = a[p, i], a[i, q]
                        a[p, i] = g - s * (h + g * tau)
                        a[i, q] = h + s * (g - h * tau)
                    for i in range(q + 1, n):
                        g, h = a[p, i], a[q, i]
                        a[p, i] = g - s * (h + g * tau)
                        a[q, i] = h + s * (g - h * tau)
                    for i in range(n):
                        g, h = v[i, p], v[i, q]
                        v[i, p] = g - s * (h + g * tau)
                        v[i, q] = h + s * (g - h * tau)

        bw += zw
        d[:] = bw
        zw[:] = 0.0

    return d.copy(), v


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
