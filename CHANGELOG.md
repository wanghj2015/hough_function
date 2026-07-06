# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/). The project is archived on
Zenodo under concept DOI
[10.5281/zenodo.21201407](https://doi.org/10.5281/zenodo.21201407).

## [1.1.0] - 2026-07-06

### Added
- Fortran finite-difference wind method (`theta_to_theta_uv.f90`) and a
  `--wind=auto|fd|groves` option. `auto` uses `fd` for `s=1`, where the Groves
  upward recurrence is numerically unstable, giving stable DW1 `u`/`v`.
- `plot_uv_modes` scripts (Python and Fortran) plotting the `U`/`V` wind modes
  for the (1,-1) and (2,2) tides, labelled `U[2,2]/3`, `U[1,-1]/10`.
- Cross-solver regression tests: cheb L2-normalization and cheb-vs-nalp wind
  amplitude agreement.
- `python/README.md` documenting the Python port, its scripts and tests.

### Changed
- Python Chebyshev solver (`cheb_hough`) now L2-normalizes each scalar mode
  (`int hough^2 dx = 1`), so its winds share the same physical amplitude as the
  ALP/Fortran solvers (an extension over the MATLAB source).
- Reproduced Fig. 2 wind panels are divided by the fixed `3, 9, 16` factors
  with divisor labels (`[2,2]/3`, ...), matching the published figure; compact
  figure sizes and thinner curves.
- Restructured docs: renamed `python/examples/` -> `python/scripts/` and
  `docs/reference.md` -> `docs/README.md`; slimmed the main README Python
  section.
- `pyproject.toml` sets `pythonpath` and the test file bootstraps `sys.path`,
  so bare `pytest` and direct execution of the test file both work.

### Unchanged
- Computed eigenvalues / equivalent depths are identical to v1.0.0.

## [1.0.0] - 2026-07-04

### Added
- Initial archived release: MATLAB, Python, and Fortran90 (CMake) solvers with
  Chebyshev-collocation and normalized-ALP methods reproducing Figs. 1-3 of
  Wang, Boyd & Akmaev (2016). Zenodo/CITATION metadata and DOI badge.
