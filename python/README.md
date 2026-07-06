# Computation of Hough Functions (Python)

Python port of the two solvers from Wang, Boyd & Akmaev (2016):

- **Chebyshev collocation** (`hough/cheb_hough.py`) — discretizes the tidal
  operator on a Chebyshev grid and solves a dense eigenproblem.
- **Normalized-ALP** (`hough/nalp_hough.py`) — expands in normalized associated
  Legendre polynomials, giving symmetric tridiagonal eigenproblems split into
  symmetric / anti-symmetric families.

The two agree to ~4 significant figures on the physical equivalent depths — the
central cross-check of the paper. See [`../docs/README.md`](../docs/README.md)
for background and the sign-convention notes.

## Install & quick start

```bash
cd python
pip install -r requirements.txt          # numpy, matplotlib

python -c "from hough import cheb_hough; print(cheb_hough.compute().h[:6])"
```

`compute(s, sigma, N)` returns a `HoughResult` with fields `x`, `lamb`,
`hough`, `hough_u`, `hough_v`, `h` (equivalent depth, km). `nalp_hough` also
offers `solve_parity(s, sigma)`, which returns the symmetric / anti-symmetric
families separately. Example tidal components: DW1 (`s=1, sigma=0.5`),
SW2 (`s=2, sigma=1.0`), TW3 (`s=3, sigma=1.5`).

## Scripts

Run from the `python/` directory as modules:

| Command | What it does |
|---------|--------------|
| `python -m scripts.plot_modes --method nalp --s 1 --sigma 0.5` | Plot the leading modes for any tide |
| `python -m scripts.plot_uv_modes` | U/V wind modes for (1,-1) and (2,2) → `../docs/uv_modes.png` |
| `python -m scripts.plot_paper_figures` | Reproduce Figs 1-3 of the paper → `../docs/fig*.png` |

`plot_uv_modes.py` uses the Chebyshev solver, L2-normalizes each mode, and
divides each mode's `U`,`V` by a fixed label factor (SW2 `(2,2)/3`,
DW1 `(1,-1)/10`); the paper figures divide the wind panels by `3, 9, 16`.

Because it differentiates the mode with the Chebyshev **spectral** operator,
the DW1 `(1,-1)` winds come out smooth through the ±30° critical latitude
(`sin φ = σ`) and at the poles. The Fortran counterpart
(`fortran/scripts/plot_uv_modes.py`) uses a finite-difference derivative
(`--wind=fd`) and shows small kinks there; the two otherwise agree.

## Tests

`tests/test_cross_check.py` guards that the two solvers stay numerically
consistent — on eigenvalues (equivalent depths) and on normalized wind
amplitude.

| Test | Checks | Why |
|------|--------|-----|
| `test_methods_agree_on_equivalent_depths` | top-7 DW1 depths match between cheb and nalp (`rtol=1e-3`) | two independent methods must agree — the paper's central claim |
| `test_gravest_dw1_depth` | gravest DW1 depth ≈ 0.694 km | absolute sanity check vs the published value |
| `test_cheb_scalar_is_l2_normalized` | every physical cheb mode has `∫ hough² dx = 1` | guards cheb's L2-normalization |
| `test_cheb_and_nalp_wind_amplitudes_agree` | gravest SW2 zonal-wind peak matches within 2% (≈ 2.82) | winds share the same physical amplitude scale |

A helper `_physical_depths` filters out the one spurious near-infinite
eigenvalue the collocation method produces. The suite does **not** test the
plotting scripts (they just emit PNGs) or the Fortran solver (cross-checked
separately — see [`../fortran/README.md`](../fortran/README.md)).

Run (from `python/`):

```bash
pytest                              # or: python3 -m pytest tests
python3 tests/test_cross_check.py   # runs the 4 and prints "cross-check tests passed"
```

`pyproject.toml` sets `pythonpath = ["."]` so `hough` imports without an
install (makes bare `pytest` work); the test file also bootstraps `sys.path`
for direct execution.

## Implementation notes

- `utils.py` helpers (`lgwt`, `pmn_polynomial_value`, `central_diff`) were
  referenced but not shipped with the original MATLAB; they are reimplemented
  here. `utils.py` also provides `jacobi_eigenvalue`, a port of Burkardt's
  Jacobi rotation eigensolver, which `nalp_hough` uses instead of
  `numpy.linalg.eigh` for the symmetric F1/F2 matrices — the paper cites this
  algorithm specifically, and it reproduces the published figures' eigenvector
  sign convention (see [`../docs/README.md`](../docs/README.md)).
- `cheb_hough` L2-normalizes each scalar mode (`∫ hough² dx = 1`), an extension
  over the MATLAB source, so its winds share nalp's physical amplitude. nalp is
  L2-normalized for free (orthonormal ALP basis + unit-norm coefficients =
  Parseval); the Chebyshev eigenvectors are only unit *Euclidean* norm on the
  grid, so they need the explicit step.
