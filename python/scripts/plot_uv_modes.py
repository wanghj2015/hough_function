"""Plot the U/V wind Hough modes for the (1,-1) and (2,2) tides.

Uses the Chebyshev collocation solver, whose spectral derivative gives winds
that stay smooth through the DW1 critical latitude (sin(lat) = sigma, i.e.
+-30 deg) and at the poles.

Conventions (as in this repo): U = zonal wind = hough_u,
V = meridional wind = hough_v.

Normalization: cheb_hough L2-normalizes each scalar mode (int hough^2 dx = 1
over [-1,1]), so the winds are on the same physical amplitude as nalp/Fortran.
Each mode's U, V are then divided by a fixed per-mode factor (~ its peak
amplitude), shown in the label: SW2 (2,2) by 3, DW1 (1,-1) by 10. Both then
sit near unit peak while the true relative size stays recoverable from the
label (DW1 is ~10/3 larger).

The Chebyshev solver returns the physical modes interleaved with spurious
collocation eigenvalues (some with near-zero or near-infinite equivalent
depth), so there is no fixed column index for a given mode. We instead pick
the gravest mode of the relevant branch, after dropping the spurious
near-infinite-depth modes:
  (2,2)  = SW2 gravest PROPAGATING mode -> largest positive depth (~ 7.86 km).
  (1,-1) = DW1 gravest TRAPPED     mode -> most negative depth   (~ -12.14 km).

Run from the ``python/`` directory:  python -m scripts.plot_uv_modes
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless: we only save a PNG
import matplotlib.pyplot as plt

from hough import cheb_hough

DOCS = os.path.join(os.path.dirname(__file__), "..", "..", "docs")

# A depth above this is a spurious near-infinite-depth collocation mode.
SPURIOUS_KM = 1.0e6

# (s, sigma, branch, divisor, colour, tag) for each mode. "branch" picks the
# gravest mode: "propagating" = largest positive depth, "trapped" = most
# negative. "divisor" is the fixed factor each L2-normalized u/v is divided by
# (~ its peak amplitude), shown in the label -- SW2 /3, DW1 /10 -- so both sit
# near unit peak while the true relative size stays recoverable from the label.
MODES = [
    dict(s=2, sigma=1.0, branch="propagating", divisor=3.0,  color="#1f4fd8", tag="2,2"),  # SW2
    dict(s=1, sigma=0.5, branch="trapped",     divisor=10.0, color="#d81f1f", tag="1,-1"), # DW1
]

def select_mode(r, branch):
    """Column index of the gravest physical mode of the requested branch."""
    h = np.where(np.abs(r.h) < SPURIOUS_KM, r.h, np.nan)  # drop spurious inf-depth
    return int(np.nanargmax(h) if branch == "propagating" else np.nanargmin(h))


def latitude(x):
    return np.degrees(np.arcsin(np.clip(x, -1.0, 1.0)))


def profiles(r, col):
    """Full-latitude profiles: U = zonal = hough_u, V = meridional = hough_v.
    cheb_hough already L2-normalizes each mode, so the winds are on the same
    physical amplitude scale as nalp/Fortran."""
    lat = latitude(r.x)
    order = np.argsort(lat)
    u = r.hough_u[:, col][order]
    v = r.hough_v[:, col][order]
    return lat[order], u, v


def scale_pair(u, v, lat, divisor, positive_lobe):
    """Divide a mode's (L2-normalized) U, V by a fixed factor and fix its sign.

    An eigenvector has no intrinsic sign (a mode and its negative are the same
    mode), and the Chebyshev eig() returns an arbitrary one -- so we pin a
    reproducible convention: sample U at 45N and flip if its sign disagrees
    with `positive_lobe` (propagating (2,2) -> positive lobe, trapped (1,-1)
    -> negative). U and V are flipped *together*, never independently: both
    are linear in the same scalar mode, so their overall sign is one shared
    gauge while their relative sign is physical.
    """
    u, v = u / divisor, v / divisor
    if (u[np.argmin(np.abs(lat - 45.0))] > 0) != positive_lobe:
        u, v = -u, -v
    return u, v


def main():
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    for m in MODES:
        r = cheb_hough.compute(s=m["s"], sigma=m["sigma"])
        col = select_mode(r, m["branch"])
        lat, u, v = profiles(r, col)
        u, v = scale_pair(u, v, lat, m["divisor"],
                          positive_lobe=(m["branch"] == "propagating"))
        d = int(m["divisor"])
        ax.plot(lat, u, color=m["color"], lw=2,
                label=rf"$U[{m['tag']}]/{d}$")
        ax.plot(lat, v, color=m["color"], lw=2, ls="--",
                label=rf"$V[{m['tag']}]/{d}$")
        print(f"({m['tag']}): cheb mode h = {r.h[col]:.4f} km")

    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlim(-90, 90)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(["90S", "60S", "30S", "0", "30N", "60N", "90N"])
    ax.set_xlabel("latitude")
    ax.set_ylabel("wind modes (normalized)")
    ax.grid(True, ls="--", alpha=0.4)
    ax.set_title(r"Hough wind modes $U,V$ for the (1,-1) and (2,2) tides",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    fig.tight_layout()

    os.makedirs(DOCS, exist_ok=True)
    path = os.path.abspath(os.path.join(DOCS, "uv_modes.png"))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)
    return path


if __name__ == "__main__":
    main()
