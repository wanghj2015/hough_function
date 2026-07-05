# Hough functions (Fortran)

Fortran90 solver for Hough functions from:

Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. *Geoscientific Model Development*, 9(4), 1477-1488.
See `../docs/reference.md` for background and the sign-convention notes.

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

### Comparing eigensolvers

`hough_coef` diagonalizes two symmetric tridiagonal matrices (F1, F2) per
run. The default is Jacobi rotations (Burkardt 2013), which is what the
paper cites and the only method guaranteed to reproduce its published
eigenvector sign convention (see `../docs/reference.md`). To sanity-check
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
