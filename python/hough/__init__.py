"""Compute Hough functions two ways, following Wang, Boyd & Akmaev (2016).

- ``cheb_hough.compute``  : Chebyshev collocation method
- ``nalp_hough.compute``  : normalized associated Legendre polynomials

Reference
---------
Wang, H., Boyd, J. P., & Akmaev, R. A. (2016). On computation of Hough
functions. Geoscientific Model Development, 9(4), 1477-1488.
https://doi.org/10.5194/gmd-9-1477-2016
"""

from . import cheb_hough, nalp_hough
from .cheb_boyd import cheb_boyd
from .utils import lgwt, pmn_polynomial_value, central_diff

__all__ = [
    "cheb_hough",
    "nalp_hough",
    "cheb_boyd",
    "lgwt",
    "pmn_polynomial_value",
    "central_diff",
]
