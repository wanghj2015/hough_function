"""Reproduce Figures 1, 2 and 3 of Wang, Boyd & Akmaev (2016).

Figures are computed with the normalized-ALP method (``nalp_hough``), as in
the paper, and written as PNGs into the repository ``docs/`` folder.

    Fig 1  DW1 (s=1, sigma=0.5)  - scalar Hough modes, 4 panels
    Fig 2  SW2 (s=2, sigma=1.0)  - scalar + u + v wind, 6 panels
    Fig 3  TW3 (s=3, sigma=1.5)  - scalar Hough modes, 2 panels

Run from the ``python/`` directory:  python -m scripts.plot_paper_figures
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless: we only save PNGs
import matplotlib.pyplot as plt

from hough import nalp_hough

# paper colour / style cycle: blue solid, red dashed, green dotted
STYLES = [
    dict(color="#1f4fd8", ls="-", lw=1.2),
    dict(color="#d81f1f", ls="--", lw=1.2),
    dict(color="#1a8a2a", ls=":", lw=1.5),
]

DOCS = os.path.join(os.path.dirname(__file__), "..", "..", "docs")


def latitude(x):
    """Convert x = sin(latitude) to latitude in degrees."""
    return np.degrees(np.arcsin(np.clip(x, -1.0, 1.0)))


def _split(family):
    """Indices of a ParityResult split into (inf, positive, negative) modes.

    Both families are returned gravest-first. Modes are sorted descending by
    equivalent depth h (== eigenvalue). The gravest positive (propagating)
    mode has the largest +h; the gravest negative (Rossby) mode has the
    largest |-h| (most negative), so the negative list is reversed.
    """
    h = family.h
    inf = list(np.where(~np.isfinite(h))[0])
    pos = list(np.where(np.isfinite(h) & (h > 0))[0])         # gravest first
    neg = list(np.where(np.isfinite(h) & (h < 0))[0])[::-1]   # most -ve first
    return inf, pos, neg


# Fixed per-mode divisors for the wind panels (compare_hough_uv.ncl
# convention): the k-th mode of a family is divided by ~its L2-normalized peak
# amplitude, so the family is shown at consistent relative scale rather than
# each curve rescaled to unit peak.
WIND_DIVISORS = [3.0, 9.0, 16.0]


def _plot_panel(ax, x, cols, labels, tag, ylim=None, normalize=False,
                divisors=None):
    lat = latitude(x)
    order = np.argsort(lat)
    lat = lat[order]
    for i, (col, lab, st) in enumerate(zip(cols, labels, STYLES)):
        # Sign is already locked in solve_parity (leading-coefficient
        # convention); do not re-flip here.
        if divisors is not None:
            y = col / divisors[i]                    # fixed family divisor
        elif normalize:
            y = col / np.max(np.abs(col))            # unit peak
        else:
            y = col
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
    sol = nalp_hough.solve_parity(s=1, sigma=0.5)
    sym, anti = sol.symmetric, sol.antisymmetric
    _, sp, sn = _split(sym)
    ai, ap, an = _split(anti)

    fig, ax = plt.subplots(2, 2, figsize=(7, 5))

    # (a) symmetric: [-1] [+1] [+3]
    _plot_panel(ax[0, 0], sym.x,
                [sym.hough[:, sn[0]], sym.hough[:, sp[0]], sym.hough[:, sp[1]]],
                ["[-1]", "[+1]", "[+3]"], "a")
    # (b) symmetric: [-1] [-3] [-5]
    _plot_panel(ax[0, 1], sym.x,
                [sym.hough[:, sn[0]], sym.hough[:, sn[1]], sym.hough[:, sn[2]]],
                ["[-1]", "[-3]", "[-5]"], "b")
    # (c) anti-symmetric: [0] [+2] [+4]
    _plot_panel(ax[1, 0], anti.x,
                [anti.hough[:, ai[0]], anti.hough[:, ap[0]], anti.hough[:, ap[1]]],
                ["[0]", "[+2]", "[+4]"], "c")
    # (d) anti-symmetric: [0] [-2] [-4]
    _plot_panel(ax[1, 1], anti.x,
                [anti.hough[:, ai[0]], anti.hough[:, an[0]], anti.hough[:, an[1]]],
                ["[0]", "[-2]", "[-4]"], "d")

    fig.suptitle(r"Figure 1.  DW1 Hough modes  ($s=1,\ \sigma=0.5$)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return _save(fig, "fig1_dw1.png")


def figure2():
    """SW2 scalar + zonal + meridional wind, 6 panels."""
    sol = nalp_hough.solve_parity(s=2, sigma=1.0)
    sym, anti = sol.symmetric, sol.antisymmetric
    # first three propagating (positive-h) modes of each family
    _, sp, _ = _split(sym)
    _, ap, _ = _split(anti)
    sp, ap = sp[:3], ap[:3]

    sd = sym.degree[sp]        # degrees 2,4,6
    ad = anti.degree[ap]       # degrees 3,5,7
    slab = [f"[2,{d}]" for d in sd]
    alab = [f"[2,{d}]" for d in ad]
    # wind labels carry the divisor, e.g. [2,2]/3 (as in the paper's Fig 2)
    slab_w = [f"[2,{d}]/{int(v)}" for d, v in zip(sd, WIND_DIVISORS)]
    alab_w = [f"[2,{d}]/{int(v)}" for d, v in zip(ad, WIND_DIVISORS)]

    fig, ax = plt.subplots(3, 2, figsize=(7, 7.5))
    # (a,b) scalar
    _plot_panel(ax[0, 0], sym.x, [sym.hough[:, i] for i in sp], slab, "a")
    _plot_panel(ax[0, 1], anti.x, [anti.hough[:, i] for i in ap], alab, "b")
    # (c,d) zonal wind u   (divided by the fixed family factors 3, 9, 16)
    _plot_panel(ax[1, 0], sym.x, [sym.hough_u[:, i] for i in sp], slab_w, "c",
                divisors=WIND_DIVISORS)
    _plot_panel(ax[1, 1], anti.x, [anti.hough_u[:, i] for i in ap], alab_w, "d",
                divisors=WIND_DIVISORS)
    # (e,f) meridional wind v  (parity flips: v of symmetric scalar is anti, etc.)
    _plot_panel(ax[2, 0], sym.x, [sym.hough_v[:, i] for i in sp], slab_w, "e",
                divisors=WIND_DIVISORS)
    _plot_panel(ax[2, 1], anti.x, [anti.hough_v[:, i] for i in ap], alab_w, "f",
                divisors=WIND_DIVISORS)

    fig.suptitle(r"Figure 2.  SW2 Hough modes  ($s=2,\ \sigma=1$)"
                 "   left: symmetric   right: anti-symmetric", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return _save(fig, "fig2_sw2.png")


def figure3():
    """TW3 scalar Hough modes, 2 panels."""
    sol = nalp_hough.solve_parity(s=3, sigma=1.5)
    sym, anti = sol.symmetric, sol.antisymmetric
    _, sp, _ = _split(sym)
    _, ap, _ = _split(anti)
    sp, ap = sp[:3], ap[:3]

    fig, ax = plt.subplots(1, 2, figsize=(7, 3))
    _plot_panel(ax[0], sym.x, [sym.hough[:, i] for i in sp],
                [f"[3,{d}]" for d in sym.degree[sp]], "a")
    _plot_panel(ax[1], anti.x, [anti.hough[:, i] for i in ap],
                [f"[3,{d}]" for d in anti.degree[ap]], "b")
    fig.suptitle(r"Figure 3.  TW3 Hough modes  ($s=3,\ \sigma=1.5$)"
                 "   left: symmetric   right: anti-symmetric", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return _save(fig, "fig3_tw3.png")


def _save(fig, name):
    os.makedirs(DOCS, exist_ok=True)
    path = os.path.abspath(os.path.join(DOCS, name))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)
    return path


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
