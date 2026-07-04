"""Compute Hough functions using normalized associated Legendre
polynomials (ALP).

Python port of nalp_hough.m from:
Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. Geoscientific Model Development, 9(4), 1477-1488.
"""

from collections import namedtuple

import numpy as np

from .utils import lgwt, pmn_polynomial_value, central_diff, jacobi_eigenvalue

# Physical constants
A_EARTH = 6.370e6                       # Earth radius (m)
G = 9.81                                # gravity (m s^-2)
OMEGA = 2.0 * np.pi / (24.0 * 3600.0)   # Earth rotation rate (s^-1)

HoughResult = namedtuple(
    "HoughResult", ["x", "lamb", "hough", "hough_u", "hough_v", "h"])


def compute(s=1, sigma=0.5, N=62, nlat=94):
    """Compute Hough functions via normalized associated Legendre polynomials.

    Parameters
    ----------
    s : int
        Zonal wavenumber. Common cases: DW1 (s=1, sigma=0.5),
        SW2 (s=2, sigma=1.0), TW3 (s=3, sigma=1.5).
    sigma : float
        Normalized frequency (frequency / (2*OMEGA)).
    N : int
        Spectral truncation (number of Legendre polynomials). Must be even.
    nlat : int
        Number of Legendre-Gauss quadrature latitudes.

    Returns
    -------
    HoughResult
        Fields: x (grid), lamb (eigenvalues), hough (mode shapes,
        columns are modes), hough_u / hough_v (wind components),
        h (equivalent depth, km).
    """
    s = int(s)
    N2 = N // 2
    # NumPy float so a zero denominator yields Inf (caught below) as in MATLAB.
    sf = np.float64(s) / np.float64(sigma)

    # define L(r) and M(r)
    L = np.zeros(N)
    M = np.zeros(N)
    with np.errstate(divide="ignore", invalid="ignore"):
        for r in range(s, N + s):
            i = r - s  # 0-based index
            L[i] = (np.sqrt((r + s + 1) * (r + s + 2)
                            * (r - s + 1) * (r - s + 2))
                    / ((2 * r + 3) * np.sqrt((2 * r + 1) * (2 * r + 5))
                       * (sf - (r + 1) * (r + 2))))
            if (s == 2) and (r == 2):
                M[i] = (-(sigma ** 2 * (sf - r * (r + 1)))
                        / ((r * (r + 1)) ** 2)
                        + (r + 2) ** 2 * (r + s + 1) * (r - s + 1)
                        / ((r + 1) ** 2 * (2 * r + 3) * (2 * r + 1)
                           * (sf - (r + 1) * (r + 2))))
            else:
                M[i] = (-(sigma ** 2 * (sf - r * (r + 1)))
                        / ((r * (r + 1)) ** 2)
                        + (r + 2) ** 2 * (r + s + 1) * (r - s + 1)
                        / ((r + 1) ** 2 * (2 * r + 3) * (2 * r + 1)
                           * (sf - (r + 1) * (r + 2)))
                        + (r - 1) ** 2 * (r ** 2 - s ** 2)
                        / (r ** 2 * (4 * r ** 2 - 1) * (sf - r * (r - 1))))
            if np.isinf(M[i]):
                M[i] = np.finfo(float).max

    # build F1 & F2 matrix
    f1 = np.zeros((N2, N2))
    f2 = np.zeros((N2, N2))
    for i in range(N2):
        f1[i, i] = M[2 * i]      # MATLAB M(2*i-1), 1-based
        f2[i, i] = M[2 * i + 1]  # MATLAB M(2*i)
        if i + 1 < N2:
            f1[i, i + 1] = L[2 * i]
            f1[i + 1, i] = L[2 * i]
            f2[i, i + 1] = L[2 * i + 1]
            f2[i + 1, i] = L[2 * i + 1]

    # symmetric modes -- solved via Jacobi rotations, as the paper specifies
    # (Burkardt 2013), not a generic LAPACK eig/eigh call: the two give the
    # same eigenvalues but a different (still valid) eigenvector sign gauge,
    # and only the Jacobi gauge reproduces the published figures. See
    # docs/reference.md.
    lamb1, v1 = jacobi_eigenvalue(f1)
    ii = np.argsort(-lamb1)
    lamb1 = lamb1[ii]
    v1 = v1[:, ii]

    # anti-symmetric modes
    lamb2, v2 = jacobi_eigenvalue(f2)
    ii = np.argsort(-lamb2)
    lamb2 = lamb2[ii]
    v2 = v2[:, ii]

    # Legendre-Gauss quadrature points
    x, w = lgwt(nlat, -1, 1)

    # normalized associated Legendre functions
    prs = pmn_polynomial_value(nlat, N + s, s, x)

    # compute Hough modes
    h1 = np.zeros((nlat, N2))
    h2 = np.zeros((nlat, N2))
    for i in range(N2):
        for j in range(N2):
            # 0-based prs column = Legendre degree. Symmetric family uses
            # degrees s, s+2, ...; anti-symmetric uses s+1, s+3, ...
            # (MATLAB used 1-based columns 2*j+s-1 / 2*j+s with j from 1.)
            i1 = 2 * j + s
            i2 = 2 * j + s + 1
            h1[:, i] += v1[j, i] * prs[:, i1]
            h2[:, i] += v2[j, i] * prs[:, i2]

    # put them together
    lamb = np.zeros(N)
    hough = np.zeros((nlat, N))
    for i in range(N2):
        lamb[2 * i] = lamb1[i]
        lamb[2 * i + 1] = lamb2[i]
        hough[:, 2 * i] = h1[:, i]
        hough[:, 2 * i + 1] = h2[:, i]

    ii = np.argsort(1.0 / lamb)
    lamb = lamb[ii]
    hough = hough[:, ii]

    # equivalent depth (km); the realmax guard above yields one spurious
    # overflow-to-inf mode, which downstream code filters out.
    with np.errstate(over="ignore"):
        h = 4.0 * A_EARTH ** 2 * OMEGA ** 2 / G * lamb / 1000.0

    # Hough functions for wind components
    b1 = (sigma ** 2 - x ** 2) * np.sqrt(1.0 - x ** 2)
    b2 = np.sqrt(1.0 - x ** 2) / (sigma ** 2 - x ** 2)
    dhdx = central_diff(hough, x)
    hough_u = (np.diag(s / b1) @ hough
               - np.diag(b2 * x / sigma) @ dhdx)
    hough_v = (np.diag((s / sigma) * x / b1) @ hough
               - np.diag(b2) @ dhdx)

    return HoughResult(x, lamb, hough, hough_u, hough_v, h)


ParityResult = namedtuple(
    "ParityResult", ["x", "lamb", "hough", "hough_u", "hough_v", "h", "degree"])
ParitySolution = namedtuple("ParitySolution", ["symmetric", "antisymmetric"])


def _winds(hough, x, s, sigma):
    """Zonal/meridional wind Hough components for a set of mode columns."""
    b1 = (sigma ** 2 - x ** 2) * np.sqrt(1.0 - x ** 2)
    b2 = np.sqrt(1.0 - x ** 2) / (sigma ** 2 - x ** 2)
    dhdx = central_diff(hough, x)
    u = np.diag(s / b1) @ hough - np.diag(b2 * x / sigma) @ dhdx
    v = np.diag((s / sigma) * x / b1) @ hough - np.diag(b2) @ dhdx
    return u, v


def solve_parity(s=1, sigma=0.5, N=62, nlat=94):
    """Solve for Hough modes keeping the two parity families separate.

    Returns a ParitySolution(symmetric, antisymmetric); each field is a
    ParityResult whose columns are modes ordered by descending eigenvalue
    (i.e. descending equivalent depth). ``degree`` holds the associated
    Legendre degree n = s + 2k (symmetric) or s + 1 + 2k (anti-symmetric)
    conventionally attached to the k-th mode of that family.

    This shares the eigenproblem construction with :func:`compute` but does
    not interleave the families, which is what the paper's figures need.
    """
    s = int(s)
    N2 = N // 2
    sf = np.float64(s) / np.float64(sigma)
    x, _ = lgwt(nlat, -1, 1)

    L = np.zeros(N)
    M = np.zeros(N)
    with np.errstate(divide="ignore", invalid="ignore"):
        for r in range(s, N + s):
            i = r - s
            L[i] = (np.sqrt((r + s + 1) * (r + s + 2)
                            * (r - s + 1) * (r - s + 2))
                    / ((2 * r + 3) * np.sqrt((2 * r + 1) * (2 * r + 5))
                       * (sf - (r + 1) * (r + 2))))
            if (s == 2) and (r == 2):
                M[i] = (-(sigma ** 2 * (sf - r * (r + 1))) / ((r * (r + 1)) ** 2)
                        + (r + 2) ** 2 * (r + s + 1) * (r - s + 1)
                        / ((r + 1) ** 2 * (2 * r + 3) * (2 * r + 1)
                           * (sf - (r + 1) * (r + 2))))
            else:
                M[i] = (-(sigma ** 2 * (sf - r * (r + 1))) / ((r * (r + 1)) ** 2)
                        + (r + 2) ** 2 * (r + s + 1) * (r - s + 1)
                        / ((r + 1) ** 2 * (2 * r + 3) * (2 * r + 1)
                           * (sf - (r + 1) * (r + 2)))
                        + (r - 1) ** 2 * (r ** 2 - s ** 2)
                        / (r ** 2 * (4 * r ** 2 - 1) * (sf - r * (r - 1))))
            if np.isinf(M[i]):
                M[i] = np.finfo(float).max

    f1 = np.zeros((N2, N2))
    f2 = np.zeros((N2, N2))
    for i in range(N2):
        f1[i, i] = M[2 * i]
        f2[i, i] = M[2 * i + 1]
        if i + 1 < N2:
            f1[i, i + 1] = f1[i + 1, i] = L[2 * i]
            f2[i, i + 1] = f2[i + 1, i] = L[2 * i + 1]

    prs = pmn_polynomial_value(nlat, N + s, s, x)

    def family(fmat, degree0):
        # degree0 = Legendre degree of the gravest mode of this family
        # (s for symmetric, s+1 for anti-symmetric).
        #
        # Solved via Jacobi rotations (Burkardt 2013), matching the method
        # the paper specifies for these symmetric matrices, rather than a
        # generic LAPACK eig/eigh call. Eigenvector sign is only defined up
        # to a gauge; eig/eigh (MATLAB and NumPy agree with each other here,
        # verified against a real MATLAB R2024a run) pick a different, only
        # sometimes-matching gauge than the published PDF. Jacobi rotations
        # -- which build the eigenvectors explicitly from the identity
        # through a deterministic sequence of plane rotations -- reproduce
        # the paper's sign for essentially every mode checked by pixel-level
        # digitization of Figs. 1-3. See docs/reference.md.
        lamb, v = jacobi_eigenvalue(fmat)
        order = np.argsort(-lamb)
        lamb, v = lamb[order], v[:, order]
        hough = np.zeros((nlat, N2))
        for i in range(N2):
            for j in range(N2):
                hough[:, i] += v[j, i] * prs[:, degree0 + 2 * j]
        with np.errstate(over="ignore"):
            h = 4.0 * A_EARTH ** 2 * OMEGA ** 2 / G * lamb / 1000.0

        u, vwind = _winds(hough, x, s, sigma)
        degree = degree0 + 2 * np.arange(N2)  # conventional degree per mode
        return ParityResult(x, lamb, hough, u, vwind, h, degree)

    symmetric = family(f1, s)          # degrees s, s+2, ...
    antisymmetric = family(f2, s + 1)  # degrees s+1, s+3, ...
    return ParitySolution(symmetric, antisymmetric)


if __name__ == "__main__":
    res = compute()
    print("first equivalent depths (km):", np.round(res.h[:6], 4))
