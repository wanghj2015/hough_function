"""Reproduce Figures 1, 2 and 3 of Wang, Boyd & Akmaev (2016) from the
*Fortran* solver's output.

This mirrors python/examples/paper_figures.py, which reproduces the same
figures from the Python nalp_hough port -- here the data instead comes from
hough_main's binary output, as a cross-check that the CLI/CMake rewrite
still reproduces the paper.

Usage:
    ./scripts/run_paper_cases.sh          # builds + runs dw1, sw2, tw3
    python3 scripts/plot_paper_figures.py # reads output/*.dat, writes PNGs
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless: we only save PNGs
import matplotlib.pyplot as plt
from scipy.io import FortranFile

HERE = os.path.dirname(os.path.abspath(__file__))
FORTRAN_DIR = os.path.dirname(HERE)
OUTPUT_DIR = os.path.join(FORTRAN_DIR, "output")

# Physical constants (must match src/hough_main.f90)
A_EARTH = 6.370e6
G = 9.81
OMEGA = 2.0 * np.pi / (24.0 * 3600.0)

# paper colour / style cycle: blue solid, red dashed, green dotted
STYLES = [
    dict(color="#1f4fd8", ls="-", lw=1.8),
    dict(color="#d81f1f", ls="--", lw=1.8),
    dict(color="#1a8a2a", ls=":", lw=2.2),
]


def latitude(x):
    """Convert x = sin(latitude) to latitude in degrees."""
    return np.degrees(np.arcsin(np.clip(x, -1.0, 1.0)))


def read_case(filename, s):
    """Read one hough_main .dat file, split into symmetric/anti-symmetric families.

    hough_main writes: nlat, nn ; latd(nlat) ; lambda(nn) ; theta(nlat,nn) ;
    theta_u(nlat,nn) ; theta_v(nlat,nn) -- lambda/theta/theta_u/theta_v pack
    the symmetric family into columns [0:n2) and the anti-symmetric family
    into columns [n2:nn), each sorted by descending eigenvalue (Jacobi
    convention -- see docs/reference.md).
    """
    path = os.path.join(OUTPUT_DIR, filename)
    with FortranFile(path, "r") as f:
        nlat, nn = f.read_ints(dtype=np.int32)
        latd = f.read_reals(dtype=np.float64)
        lamb = f.read_reals(dtype=np.float64)
        theta = f.read_reals(dtype=np.float64).reshape((nlat, nn), order="F")
        theta_u = f.read_reals(dtype=np.float64).reshape((nlat, nn), order="F")
        theta_v = f.read_reals(dtype=np.float64).reshape((nlat, nn), order="F")

    n2 = nn // 2
    x = np.sin(np.radians(latd))

    def family(idx, degree0):
        lm = lamb[idx]
        with np.errstate(over="ignore", invalid="ignore"):
            h = 4.0 * A_EARTH ** 2 * OMEGA ** 2 / G * lm / 1000.0
        return dict(
            x=x,
            lamb=lm,
            h=h,
            hough=theta[:, idx],
            hough_u=theta_u[:, idx],
            hough_v=theta_v[:, idx],
            degree=degree0 + 2 * np.arange(n2),
        )

    sym = family(slice(0, n2), s)          # degrees s, s+2, ...
    anti = family(slice(n2, nn), s + 1)    # degrees s+1, s+3, ...
    return sym, anti


def _split(family):
    """Indices of a family split into (inf, positive, negative) modes.

    Both branches are gravest-first. The gravest positive (propagating)
    mode has the largest +h; the gravest negative (Rossby) mode has the
    largest |-h| (most negative), so the negative list is reversed.
    """
    h = family["h"]
    inf = list(np.where(~np.isfinite(h))[0])
    pos = list(np.where(np.isfinite(h) & (h > 0))[0])         # gravest first
    neg = list(np.where(np.isfinite(h) & (h < 0))[0])[::-1]   # most -ve first
    return inf, pos, neg


def _plot_panel(ax, x, cols, labels, tag, ylim=None, normalize=False):
    lat = latitude(x)
    order = np.argsort(lat)
    lat = lat[order]
    for col, lab, st in zip(cols, labels, STYLES):
        y = col / np.max(np.abs(col)) if normalize else col   # unit peak (winds)
        ax.plot(lat, y[order], label=lab, **st)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(["90S", "60S", "30S", "0", "30N", "60N", "90N"],
                       fontsize=8)
    if ylim:
        ax.set_ylim(*ylim)
    ax.grid(True, ls="--", alpha=0.4)
    ax.text(0.03, 0.90, f"({tag})", transform=ax.transAxes,
            fontsize=11, fontweight="bold",
            bbox=dict(boxstyle="square", fc="white", ec="0.5"))
    ax.legend(loc="lower right", fontsize=7, framealpha=0.9)


# --------------------------------------------------------------------------
def figure1():
    """DW1 scalar Hough modes: signed-index labelling, 4 panels."""
    sym, anti = read_case("hough_s1_sigma.5000.dat", s=1)
    _, sp, sn = _split(sym)
    ai, ap, an = _split(anti)

    fig, ax = plt.subplots(2, 2, figsize=(10, 7))

    # (a) symmetric: [-1] [+1] [+3]
    _plot_panel(ax[0, 0], sym["x"],
                [sym["hough"][:, sn[0]], sym["hough"][:, sp[0]], sym["hough"][:, sp[1]]],
                ["[-1]", "[+1]", "[+3]"], "a")
    # (b) symmetric: [-1] [-3] [-5]
    _plot_panel(ax[0, 1], sym["x"],
                [sym["hough"][:, sn[0]], sym["hough"][:, sn[1]], sym["hough"][:, sn[2]]],
                ["[-1]", "[-3]", "[-5]"], "b")
    # (c) anti-symmetric: [0] [+2] [+4]
    _plot_panel(ax[1, 0], anti["x"],
                [anti["hough"][:, ai[0]], anti["hough"][:, ap[0]], anti["hough"][:, ap[1]]],
                ["[0]", "[+2]", "[+4]"], "c")
    # (d) anti-symmetric: [0] [-2] [-4]
    _plot_panel(ax[1, 1], anti["x"],
                [anti["hough"][:, ai[0]], anti["hough"][:, an[0]], anti["hough"][:, an[1]]],
                ["[0]", "[-2]", "[-4]"], "d")

    fig.suptitle(r"Figure 1.  DW1 Hough modes  ($s=1,\ \sigma=0.5$)  -- Fortran solver")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return _save(fig, "fig1_dw1.png")


def figure2():
    """SW2 scalar + zonal + meridional wind, 6 panels."""
    sym, anti = read_case("hough_s2_sigma1.0000.dat", s=2)
    _, sp, _ = _split(sym)
    _, ap, _ = _split(anti)
    sp, ap = sp[:3], ap[:3]

    sd = sym["degree"][sp]        # degrees 2,4,6
    ad = anti["degree"][ap]       # degrees 3,5,7
    slab = [f"[2,{d}]" for d in sd]
    alab = [f"[2,{d}]" for d in ad]

    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    # (a,b) scalar
    _plot_panel(ax[0, 0], sym["x"], [sym["hough"][:, i] for i in sp], slab, "a")
    _plot_panel(ax[0, 1], anti["x"], [anti["hough"][:, i] for i in ap], alab, "b")
    # (c,d) zonal wind u   (normalized to unit peak, as in the paper)
    _plot_panel(ax[1, 0], sym["x"], [sym["hough_u"][:, i] for i in sp], slab, "c",
                normalize=True)
    _plot_panel(ax[1, 1], anti["x"], [anti["hough_u"][:, i] for i in ap], alab, "d",
                normalize=True)
    # (e,f) meridional wind v  (parity flips: v of symmetric scalar is anti, etc.)
    _plot_panel(ax[2, 0], sym["x"], [sym["hough_v"][:, i] for i in sp], slab, "e",
                normalize=True)
    _plot_panel(ax[2, 1], anti["x"], [anti["hough_v"][:, i] for i in ap], alab, "f",
                normalize=True)

    for a, t in zip(ax[:, 0], ["scalar", "zonal wind u", "merid. wind v"]):
        a.set_ylabel(t, fontsize=9)
    fig.suptitle(r"Figure 2.  SW2 Hough modes  ($s=2,\ \sigma=1$)  -- Fortran solver"
                 "   left: symmetric   right: anti-symmetric")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return _save(fig, "fig2_sw2.png")


def figure3():
    """TW3 scalar Hough modes, 2 panels."""
    sym, anti = read_case("hough_s3_sigma1.5000.dat", s=3)
    _, sp, _ = _split(sym)
    _, ap, _ = _split(anti)
    sp, ap = sp[:3], ap[:3]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    _plot_panel(ax[0], sym["x"], [sym["hough"][:, i] for i in sp],
                [f"[3,{d}]" for d in sym["degree"][sp]], "a")
    _plot_panel(ax[1], anti["x"], [anti["hough"][:, i] for i in ap],
                [f"[3,{d}]" for d in anti["degree"][ap]], "b")
    fig.suptitle(r"Figure 3.  TW3 Hough modes  ($s=3,\ \sigma=1.5$)  -- Fortran solver"
                 "   left: symmetric   right: anti-symmetric")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return _save(fig, "fig3_tw3.png")


def _save(fig, name):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = os.path.join(OUTPUT_DIR, name)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)
    return path


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
