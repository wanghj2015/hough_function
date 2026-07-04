# Reference

Wang, H., Boyd, J. P., & Akmaev, R. A. (2016).
**On computation of Hough functions.**
*Geoscientific Model Development*, 9(4), 1477–1488.
https://doi.org/10.5194/gmd-9-1477-2016

## What the code computes

Hough functions are the eigenfunctions of Laplace's tidal equation and form
the meridional basis for atmospheric tides and planetary-scale waves. The
repository provides two independent solvers:

| Method | File(s) | Idea |
|--------|---------|------|
| Chebyshev collocation | `cheb_hough`, `cheb_boyd` | Discretize the tidal operator on a Chebyshev grid and solve a dense eigenproblem. |
| Normalized ALPs | `nalp_hough`, `utils` | Expand in normalized associated Legendre polynomials, giving a symmetric tridiagonal eigenproblem split into symmetric / anti-symmetric modes. |

Both return, per mode: the eigenvalue `lamb`, the Hough function `hough`,
the wind components `hough_u` / `hough_v`, and the equivalent depth `h` (km).

## Example tidal components

| Name | `s` | `sigma` |
|------|-----|---------|
| DW1  | 1   | 0.5     |
| SW2  | 2   | 1.0     |
| TW3  | 3   | 1.5     |

`sigma` is the frequency normalized by twice the Earth rotation rate.

## Validation

The two methods agree to ~4 significant figures on the physical equivalent
depths (e.g. the gravest DW1 mode ≈ 0.694 km). The Chebyshev solver also
returns one spurious near-infinite eigenvalue (an artifact of the
collocation discretization), which downstream code should ignore. See
`python/tests/test_cross_check.py`.

## Reproduced figures

`python/examples/paper_figures.py` regenerates Figures 1–3 of the paper with
the normalized-ALP method and writes them here:

| File | Paper figure | Content |
|------|--------------|---------|
| `fig1_dw1.png` | Fig. 1 | DW1 (s=1, σ=0.5) scalar modes, 4 panels, signed-index labels `[±n]`, `[0]` = missing mode |
| `fig2_sw2.png` | Fig. 2 | SW2 (s=2, σ=1) scalar + zonal `u` + meridional `v`, 6 panels, `(s,n)` labels |
| `fig3_tw3.png` | Fig. 3 | TW3 (s=3, σ=1.5) scalar modes, 2 panels |

Notes:
- Mode labelling: positive index = propagating (positive equivalent depth,
  gravest first), negative index = Rossby (most negative depth first), `[0]` =
  the infinite-depth "missing" mode; odd indices are symmetric, even
  anti-symmetric.
- Wind panels are normalized to unit peak (as in the paper); scalar panels use
  the natural L² normalization.
- **Sign convention.** An eigenvector is only defined up to a sign, and that
  sign is a physically meaningless gauge (a Hough function and its negative are
  the same mode). Both `compute()` and `solve_parity()` solve the symmetric
  F1/F2 matrices with **Jacobi rotations** (`hough.utils.jacobi_eigenvalue`, a
  faithful port of Burkardt's `jacobi_eigenvalue.m`), not `numpy.linalg.eigh`
  — matching the algorithm the paper's text actually cites for these matrices
  (Sect. 2.1: "The Fortran 90 source code of the Jacobi eigenvalue algorithm
  implemented by Burkardt (2013) can be used to solve the two symmetric
  matrix eigenvalue problems"), even though the MATLAB listing in Appendix B1
  shows a plain `eig()` call.

  This matters because eig/eigh and Jacobi rotations, while mathematically
  equivalent (same eigenvalues), pick *different* eigenvector sign gauges:
  eig/eigh normalize post hoc (e.g. largest-magnitude component positive),
  while Jacobi rotations build each eigenvector explicitly, sweep by sweep,
  starting from the identity matrix. We verified:
  - A real MATLAB R2024a run of the unmodified `nalp_hough.m`/`cheb_hough.m`
    matches NumPy's `eig`/`eigh` bit-for-bit (grid, eigenvalues, mode shapes,
    wind components) for DW1/SW2/TW3 — so eig/eigh's gauge is *not* a
    MATLAB-vs-NumPy platform artifact; a prior version of this note blamed a
    mismatch on Octave specifically, but real MATLAB agrees with NumPy.
  - That shared eig/eigh gauge nonetheless disagrees with the published PDF
    for roughly half the modes in Figs. 1–3 (verified by extracting each
    curve's actual pixels from the PDF by color and comparing sign at a dozen
    reference latitudes).
  - Switching to Jacobi rotations reproduces the published sign for every mode
    checked across DW1, SW2 and TW3 — confirmed against the PDF and against
    the regenerated `docs/fig*.png`.

  Earlier versions of this code (and of this note) tried to patch the
  eig/eigh output after the fact — first a "south-pole sign" heuristic, later
  a hardcoded per-mode flip table — because no simple formula was found. Both
  were replaced once it became clear the discrepancy was a solver-*algorithm*
  choice, not a mode-by-mode gauge to guess at.
- The parity solver used by the figures is
  `hough.nalp_hough.solve_parity(s, sigma, N, nlat)`, which returns the
  symmetric and anti-symmetric families separately.
