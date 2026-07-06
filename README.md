# Computation of Hough Functions

Code for computing Hough functions from the paper:

> Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
> functions. *Geoscientific Model Development*, 9(4), 1477–1488.
> https://doi.org/10.5194/gmd-9-1477-2016

This code is archived on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21201407.svg)](https://doi.org/10.5281/zenodo.21201407)

Three independent implementations are provided: the original **MATLAB**,
a **Python** port, and a **Fortran90** solver (CMake build). MATLAB/Python
each offer two solvers -- a **Chebyshev collocation** method and a
**normalized associated Legendre polynomial (ALP)** method -- which agree
to ~4 significant figures on the physical equivalent depths; the Fortran
solver uses the normalized-ALP approach and reproduces the paper's
published figures directly.

## Layout

```
matlab/                 Original MATLAB implementation
  cheb_boyd.m           Chebyshev differentiation matrices
  cheb_hough.m          Chebyshev collocation solver
  nalp_hough.m          Normalized-ALP solver
python/                 Python port (see python/README.md)
  hough/                Importable package
    cheb_boyd.py        Chebyshev differentiation matrices
    cheb_hough.py       Chebyshev collocation solver  (compute())
    nalp_hough.py       Normalized-ALP solver         (compute())
    utils.py            lgwt, pmn_polynomial_value, central_diff
  scripts/plot_modes.py       Plot the leading modes (any tide)
  scripts/plot_uv_modes.py    Plot U/V wind modes for (1,-1) and (2,2)
  scripts/plot_paper_figures.py  Reproduce Figs 1-3 of the paper into docs/
  tests/test_cross_check.py    Assert both methods agree
fortran/                Fortran90 solver (CMake build)
  CMakeLists.txt        Build config; auto-detects LAPACK
  src/hough_main.f90    CLI driver (--preset, --solver, --compare-solvers, ...)
  src/eigensolvers.f90  Jacobi (default) + LAPACK dstev/dsyev/dsyevd, cross-checked
  scripts/run_paper_cases.sh      Build + run dw1/sw2/tw3
  scripts/plot_paper_figures.py   Reproduce Figs 1-3 from the Fortran output
docs/README.md       Background and parameter reference
docs/fig1_dw1.png       Reproduced paper figures
docs/fig2_sw2.png
docs/fig3_tw3.png
docs/gmd_9_1477_2016.pdf   The paper
```

## Python usage

```bash
cd python
pip install -r requirements.txt

python -c "from hough import cheb_hough; print(cheb_hough.compute().h[:6])"
python -m scripts.plot_uv_modes         # U/V wind modes into ../docs/uv_modes.png
python -m scripts.plot_paper_figures    # write Figs 1-3 into ../docs/
pytest                                   # cross-check the two solvers
```

`plot_uv_modes` uses the Chebyshev **spectral** derivative, so the DW1 `(1,-1)`
winds stay smooth through the ±30° critical latitude and the poles — smoother
than the Fortran `--wind=fd` finite-difference version, which kinks slightly
there.

See [`python/README.md`](python/README.md) for the full script and test
reference, and the implementation notes (`utils.py` helpers, the Jacobi
eigensolver, and cheb's L2-normalization).

## Fortran usage

```bash
cd fortran
cmake -B build -S . && cmake --build build

./build/hough_main --preset=dw1          # or sw2, tw3
./build/hough_main -s 2 -f 1.0 --solver=dstev --compare-solvers
./build/hough_main --help

./scripts/run_paper_cases.sh             # build + run dw1/sw2/tw3
python3 scripts/plot_paper_figures.py    # reproduce Figs 1-3 into output/
```

See [`fortran/README.md`](fortran/README.md) for details, including the
`--compare-solvers` cross-check against LAPACK's dstev/dsyev/dsyevd.

## Licence

CC BY 4.0 — see [`LICENSE`](LICENSE) and the publisher's
[copyright & licence policy](https://www.geoscientific-model-development.net/policies/licence_and_copyright.html).
