"""Plot the leading Hough modes computed by either method.

Usage:
    python -m examples.plot_modes            # Chebyshev method, DW1
    python -m examples.plot_modes --method nalp --s 2 --sigma 1.0

Run from the ``python/`` directory (so ``hough`` is importable).
"""

import argparse

import matplotlib.pyplot as plt

from hough import cheb_hough, nalp_hough


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=["cheb", "nalp"], default="cheb")
    parser.add_argument("--s", type=float, default=1.0, help="zonal wavenumber")
    parser.add_argument("--sigma", type=float, default=0.5,
                        help="normalized frequency")
    parser.add_argument("--nmodes", type=int, default=60,
                        help="number of modes to plot")
    args = parser.parse_args()

    module = cheb_hough if args.method == "cheb" else nalp_hough
    res = module.compute(s=args.s, sigma=args.sigma)

    nmodes = min(args.nmodes, res.hough.shape[1])
    ncols = 6
    nrows = (nmodes + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3 * nrows))
    for j, ax in enumerate(axes.ravel()):
        if j < nmodes:
            ax.plot(res.x, res.hough[:, j], linewidth=2)
            ax.grid(True)
        else:
            ax.axis("off")
    fig.suptitle(f"{args.method} Hough modes  (s={args.s}, sigma={args.sigma})")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
