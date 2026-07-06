# Computation of Hough Functions (Fortran)

Fortran90 solver for Hough functions from:

Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. *Geoscientific Model Development*, 9(4), 1477-1488.
See `../docs/README.md` for background and the sign-convention notes.

## Build

```sh
cmake -B build -S .
cmake --build build
```

LAPACK is detected automatically (`find_package(LAPACK)`). If found, three
extra eigensolvers (`dstev`, `dsyev`, `dsyevd`) are enabled alongside the
default Jacobi-rotation solver; if not, the build still works with Jacobi
only.

## Run

```sh
./build/hough_main --preset=dw1      # s=1, sigma=0.5
./build/hough_main --preset=sw2      # s=2, sigma=1.0
./build/hough_main --preset=tw3      # s=3, sigma=1.5
./build/hough_main -s 2 -f 1.0       # or set s / sigma directly
./build/hough_main --help
```

Each run prints a short summary (parameters, normalization checks, gravest
equivalent depths) and writes an unformatted binary file to `output/`
containing the latitudes, eigenvalues and scalar/wind Hough modes.

### Wind (u,v) method: `--wind`

The u/v modes can be computed two ways:

- **`groves`** — the Groves (1981) coefficient recurrence (`hough_mode_uv`).
  Fast, but the *upward* recurrence is numerically unstable for `s=1`: it
  amplifies the high-degree eigenvector tail into a grid-scale sawtooth
  (roughness grows with `--nn`), so DW1 winds come out corrupted.
- **`fd`** — differentiate the (well-resolved) scalar mode directly
  (`theta_to_theta_uv`). Stable for all `s`. Validated to reproduce the
  legacy `theta_to_theta_uv` output to ~1e-7 and to agree with `groves` to
  ~1% for SW2 (where groves is reliable).

```sh
./build/hough_main --preset=dw1 --wind=fd       # stable DW1 winds
./build/hough_main --preset=sw2 --wind=groves   # force the recurrence
```

`--wind=auto` (the default) picks `fd` for `s=1` and `groves` otherwise, so
the scalar modes and the SW2/TW3 paper figures are unchanged.

### Comparing eigensolvers

`hough_coef` diagonalizes two symmetric tridiagonal matrices (F1, F2) per
run. The default is Jacobi rotations (Burkardt 2013), which is what the
paper cites and the only method guaranteed to reproduce its published
eigenvector sign convention (see `../docs/README.md`). To sanity-check
that other solvers agree on the eigenvalues:

```sh
./build/hough_main --preset=sw2 --solver=dstev --compare-solvers
```

`--solver` picks which result is actually used downstream; `--compare-solvers`
additionally prints every available solver's max eigenvalue deviation from
Jacobi. They should all agree to machine precision -- eigenvalues are
gauge-invariant even though eigenvector signs are not. One exception: some
(s, sigma) pairs (e.g. DW1's anti-symmetric family) have a genuine infinite
equivalent-depth mode, which LAPACK's routines reject as non-finite input;
the comparison table skips them for that matrix and says so.

## Reproducing the paper's figures

```sh
./scripts/run_paper_cases.sh        # builds, then runs dw1/sw2/tw3
python3 scripts/plot_paper_figures.py
```

This writes `fig1_dw1.png`, `fig2_sw2.png`, `fig3_tw3.png` to `output/`,
reproduced entirely from this Fortran solver's output (verified to match
`../docs/fig*.png`, which are reproduced independently by the Python port).
Requires `numpy`, `scipy` and `matplotlib`.

## Plotting the U/V wind modes

```sh
python3 scripts/plot_uv_modes.py        # builds the .dat files itself
```

Writes `output/uv_modes.png`: the zonal (`U`) and meridional (`V`) wind Hough
modes for the diurnal `(1,-1)` and semi-diurnal `(2,2)` tides, on one axis.
Each mode's `U`,`V` are divided by a fixed factor shown in the label — `(2,2)`
by 3, `(1,-1)` by 10 (≈ their peak amplitudes, so both sit near unit peak
while the true ~10/3 size ratio stays readable). DW1 uses the `fd` wind method
(`--wind=auto` picks it for `s=1`, where the Groves recurrence is unstable).

This mirrors `python/scripts/plot_uv_modes.py` (which uses the Chebyshev
solver); the two agree, with the Fortran `fd` winds showing slight kinks at
the DW1 critical latitude (±30°) that the Python spectral derivative smooths.
