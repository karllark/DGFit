import argparse
import importlib.resources as importlib_resources

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm

from dgfit.obsdata import ObsData
from dgfit.dustgrains import DustGrains


def main():
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wave", default=0.1, type=float, help="lamda in A(lambda)/A(V)"
    )
    parser.add_argument(
        "-c",
        "--composition",
        choices=[
            "astro-silicates",
            "astro-carbonaceous",
            "astro-graphite",
            "PAH-Z04",
            "Graphite-Z04",
            "Silicates-Z04",
            "ACH2-Z04",
            "Silicates1-Z04",
            "Silicates2-Z04",
            "Carbonaceous-HD23",
            "AstroDust-HD23",
            "a-C-Themis",
            "a-C:H-Themis",
            "aSil-2-Themis",
        ],
        default="astro-silicates",
        help="Grain composition",
    )
    parser.add_argument(
        "--obsdata",
        type=str,
        default="none",
        help="transform to observed data grids, with the name of the observed data file as input",
    )
    parser.add_argument(
        "--everynth", type=int, default=1, help="Use every nth grain size"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    DG = DustGrains()
    ref = importlib_resources.files("dgfit") / "data"
    with importlib_resources.as_file(ref) as data_path:
        DG.from_files(
            args.composition,
            path=str(data_path) + "/indiv_grain/",
            every_nth=args.everynth,
        )

    if args.obsdata != "none":
        OD = ObsData(args.obsdata)
        new_DG = DustGrains()
        new_DG.from_object(DG, OD)
        DG = new_DG

    plot(DG, args.wave, args.composition, args.pdf, args.png, args.eps)


def plot(DG, wave, composition, pdf=False, png=False, eps=False):

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    ws_indxs = np.argsort(DG.wavelengths)
    waves = DG.wavelengths[ws_indxs]

    fig, ax = plt.subplots(figsize=(8, 6))

    num_segments = DG.n_sizes
    DG.sizes *= 10**4
    cmap = get_cmap("hsv", num_segments)
    norm = LogNorm(vmin=min(DG.sizes), vmax=max(DG.sizes))
    colors = [cmap(i) for i in range(num_segments)]
    for i in range(DG.n_sizes):
        pcolor = colors[i]

        # get the values at specified lambda and V
        al = np.interp([wave, 0.55, 0.45], waves, DG.cext[i, ws_indxs])

        rv = al[1] / (al[2] - al[1])
        ax.plot(rv, al[0] / al[1], "o", color=pcolor)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.05, pad=0.04)
    cbar.set_label(r"Grain sizes [$\mu m$]")

    ax.set_xlabel(r"R(V)")
    ax.set_xlim(0.0, 10.0)
    ax.set_ylabel(f"A({wave})/A(V)")
    ax.set_yscale("log")
    ax.set_title(composition)

    plt.tight_layout()

    # show or save
    basename = "DustGrains_diag_%s" % (composition)
    if png:
        plt.savefig(basename + ".png")
    elif eps:
        plt.savefig(basename + ".eps")
    elif pdf:
        plt.savefig(basename + ".pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()
