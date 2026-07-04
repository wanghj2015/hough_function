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
  the same mode). `solve_parity` fixes each mode's outermost lobe (its approach
  to the south pole) to the natural pole sign of that family's leading
  associated Legendre function, `sign(P^s_{degree0}) = (-1)^degree0`. This is a
  single deterministic rule applied uniformly — not per-curve flipping.

  What it does and does not match in the paper:
  - **Symmetric** modes: matches the published signs in every panel checked
    (Figs 1a/1b, 2a, 3a) — their pole signs really do follow `(-1)^s`.
  - **Anti-symmetric** modes: matches the low-order / propagating modes, but
    some higher modes differ from the PDF (e.g. `[+2]`,`[+4]` in Fig 1c;
    `[2,7]` in Fig 2b/d/f). The paper's anti-symmetric signs follow *no*
    consistent rule (DW1 `+,-,-`; SW2 `-,-,+`; TW3 `+,+,+`) — they are just
    MATLAB `eig`'s raw output — so no clean convention can reproduce them all.

  Shapes, node counts and amplitudes are correct throughout; only some
  anti-symmetric curve *signs* differ from the PDF, which is immaterial.
  (`compute()` imposes no convention at all, mirroring the MATLAB code; only
  `solve_parity`, used for the figures, applies this one.)

  The paper's signs are not reproducible by any tool but the authors' exact
  MATLAB build: running the original `nalp_hough.m` through Octave (with `eig`
  and the original indexing) yields yet another gauge — it flips a *different*
  set of curves than NumPy (e.g. Octave matches `[+3]`/`[+4]` but flips
  `[-5]`/`[+2]` in Fig 1). Eigenvector signs are LAPACK-implementation
  dependent (MATLAB != Octave != NumPy), so per-curve sign matching is
  impossible in general; the convention above is the closest deterministic,
  principled choice.
- The parity solver used by the figures is
  `hough.nalp_hough.solve_parity(s, sigma, N, nlat)`, which returns the
  symmetric and anti-symmetric families separately.
