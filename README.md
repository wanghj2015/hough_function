# hough_function

Code for computing Hough functions from the paper:

> Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
> functions. *Geoscientific Model Development*, 9(4), 1477–1488.
> https://doi.org/10.5194/gmd-9-1477-2016

Two independent solvers are provided in both MATLAB and Python: a **Chebyshev
collocation** method and a **normalized associated Legendre polynomial (ALP)**
method. They agree to ~4 significant figures on the physical equivalent depths.

## Layout

```
matlab/                 Original MATLAB implementation
  cheb_boyd.m           Chebyshev differentiation matrices
  cheb_hough.m          Chebyshev collocation solver
  nalp_hough.m          Normalized-ALP solver
python/                 Python port
  hough/                Importable package
    cheb_boyd.py        Chebyshev differentiation matrices
    cheb_hough.py       Chebyshev collocation solver  (compute())
    nalp_hough.py       Normalized-ALP solver         (compute())
    utils.py            lgwt, pmn_polynomial_value, central_diff
  examples/plot_modes.py     Plot the leading modes (any tide)
  examples/paper_figures.py  Reproduce Figs 1-3 of the paper into docs/
  tests/test_cross_check.py  Assert both methods agree
docs/reference.md       Background and parameter reference
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
python -m examples.plot_modes --method nalp --s 1 --sigma 0.5
python -m examples.paper_figures         # write Figs 1-3 into ../docs/
pytest                                   # cross-check the two methods
```

`compute(s, sigma, N)` returns a `HoughResult` with fields `x`, `lamb`,
`hough`, `hough_u`, `hough_v`, `h` (equivalent depth in km). Example tidal
components: DW1 (`s=1, sigma=0.5`), SW2 (`s=2, sigma=1.0`), TW3 (`s=3, sigma=1.5`).

> Note: the Python `utils.py` helpers (`lgwt`, `pmn_polynomial_value`,
> `central_diff`) were referenced but not shipped with the original MATLAB;
> they are reimplemented here. `utils.py` also provides `jacobi_eigenvalue`,
> a port of Burkardt's Jacobi rotation eigensolver, which `nalp_hough` uses
> instead of `numpy.linalg.eigh` for the symmetric F1/F2 matrices — the
> paper cites this algorithm specifically, and it's what reproduces the
> published figures' eigenvector sign convention (see `docs/reference.md`).

## Licence

CC BY 4.0 — see [`LICENSE`](LICENSE) and the publisher's
[copyright & licence policy](https://www.geoscientific-model-development.net/policies/licence_and_copyright.html).
