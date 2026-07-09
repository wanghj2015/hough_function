"""Plot the U/V wind Hough modes for the (1,-1) and (2,2) tides, from the
FORTRAN solver -- the counterpart of python/examples/plot_uv_modes.py.

Conventions (as in this repo): U = zonal wind = hough_u,
V = meridional wind = hough_v.

The ALP-based Fortran solver is L2-normalized by construction, so the winds
are already on the physical amplitude scale (no extra L2 step, unlike the
Chebyshev solver). Each mode's U, V are divided by a fixed per-mode factor
(~ its peak), shown in the label: SW2 (2,2) by 3, DW1 (1,-1) by 10. DW1 winds
use the fd method (auto for s=1), where the Groves recurrence is unstable.

Modes are selected by branch: (2,2) = gravest propagating (largest positive
depth), (1,-1) = gravest trapped (most negative depth).

Usage:  python3 scripts/plot_uv_modes.py   (builds the .dat files itself)
"""
import os
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_paper_figures import read_case, _split, latitude, OUTPUT_DIR, FORTRAN_DIR

HOUGH_MAIN = os.path.join(FORTRAN_DIR, "build", "hough_main")

# (s, sigma, branch, divisor, colour, tag, datfile) per mode
MODES = [
    dict(s=2, sigma=1.0, branch="pos", divisor=3.0,  color="#1f4fd8", tag="2,2",  dat="uv_sw2.dat"),
    dict(s=1, sigma=0.5, branch="neg", divisor=10.0, color="#d81f1f", tag="1,-1", dat="uv_dw1.dat"),
]


def run(s, sigma, name):
    """Run hough_main (--wind=auto: fd for s=1, groves otherwise)."""
    subprocess.run([HOUGH_MAIN, "-s", str(s), "-f", str(sigma),
                    "-o", os.path.join(OUTPUT_DIR, name)],
                   check=True, stdout=subprocess.DEVNULL)


def select(sym, branch):
    _, pos, neg = _split(sym)
    return pos[0] if branch == "pos" else neg[0]


def scale_pair(u, v, divisor):
    """Divide the (L2-normalized) U, V by a fixed factor.

    An eigenvector has no intrinsic sign, and we keep whatever sign the solver
    returns rather than forcing a convention. U and V share one scalar mode, so
    their relative sign is physical and preserved either way.
    """
    return u / divisor, v / divisor


def main():
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    for m in MODES:
        run(m["s"], m["sigma"], m["dat"])
        sym, _ = read_case(m["dat"], s=m["s"])
        col = select(sym, m["branch"])
        lat = latitude(sym["x"])
        order = np.argsort(lat)
        lat = lat[order]
        u = sym["hough_u"][:, col][order]
        v = sym["hough_v"][:, col][order]
        u, v = scale_pair(u, v, m["divisor"])
        d = int(m["divisor"])
        ax.plot(lat, u, color=m["color"], lw=2, label=rf"$U[{m['tag']}]/{d}$")
        ax.plot(lat, v, color=m["color"], lw=2, ls="--",
                label=rf"$V[{m['tag']}]/{d}$")
        print(f"({m['tag']}): h = {sym['h'][col]:.4f} km")

    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlim(-90, 90)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(["90S", "60S", "30S", "0", "30N", "60N", "90N"])
    ax.set_xlabel("latitude")
    ax.set_ylabel("wind modes (normalized)")
    ax.grid(True, ls="--", alpha=0.4)
    ax.set_title(r"Hough wind modes $U,V$ for the (1,-1) and (2,2) tides "
                 "-- Fortran solver", fontsize=9)
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    fig.tight_layout()

    path = os.path.join(OUTPUT_DIR, "uv_modes.png")
    fig.savefig(path, dpi=150)
    print("wrote", path)


if __name__ == "__main__":
    main()
