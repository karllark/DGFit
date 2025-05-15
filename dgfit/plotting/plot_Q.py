import argparse
import importlib.resources as importlib_resources

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math

from dgfit.obsdata import ObsData
from dgfit.dustgrains import DustGrains


def main():
    # commandline parser
    parser = argparse.ArgumentParser()
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
        "--everynth", type=int, default=5, help="Use every nth grain size"
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

    plot(DG, args.composition, args.png, args.eps, args.pdf)


def plot(DG, composition, png=False, eps=False, pdf=False):
    # setup the plots
    fontsize = 12
    font = {"size": fontsize}

    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(15, 10))

    DG.sizes *= 10**4  # convert to micron
    ws_indxs = np.argsort(DG.wavelengths)
    waves = DG.wavelengths[ws_indxs]

    Qsca = DG.csca[:, ws_indxs] / (math.pi * (DG.sizes[:, np.newaxis]) ** 2)
    Qabs = DG.cabs[:, ws_indxs] / (math.pi * (DG.sizes[:, np.newaxis]) ** 2)

    total_csca = np.sum(Qsca[:, ws_indxs], axis=0)
    total_cabs = np.sum(Qabs[:, ws_indxs], axis=0)

    ax[0].plot(waves, total_csca)
    ax[0].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[0].set_ylabel("Q(sca)")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")

    ax[1].plot(waves, total_cabs)
    ax[1].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[1].set_ylabel("Q(abs)")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")

    fig.suptitle(composition, fontsize=20)
    fig.subplots_adjust(top=0.92)

    # show or save
    basename = "DustGrains_diag_%s" % (composition)
    if png:
        fig.savefig(basename + ".png")
    elif eps:
        fig.savefig(basename + ".eps")
    elif pdf:
        fig.savefig(basename + ".pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()
