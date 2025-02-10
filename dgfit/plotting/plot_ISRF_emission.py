import argparse
import importlib.resources as importlib_resources

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.cm import ScalarMappable
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm

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
            "astro-PAH-ionized",
            "astro-PAH-neutral",
        ],
        default="astro-silicates",
        help="Grain composition",
    )
    parser.add_argument("--wave", default=0.1, type=float)
    parser.add_argument("--grainsize", default=0, type=int)
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

    ws_indxs = np.argsort(DG.wavelengths_emission)
    waves = DG.wavelengths_emission[ws_indxs]

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(15, 10))

    alpha = 0.25
    interpolated = (1 - alpha) * DG.emission_1 + (alpha * DG.emission_5)

    # or i in range (DG.n_sizes):
    # em_1 = np.interp(args.wave, waves, DG.emission_1[i, ws_indxs])
    # ax.plot(DG.sizes[i] * 1e4, em_1, 'o', color="blue", label="1")

    # em_2 = np.interp(args.wave, waves, DG.emission_2[i, ws_indxs])
    # ax.plot(DG.sizes[i] * 1e4, em_2, 'o', color="red", label="2")

    # em_3 = np.interp(args.wave, waves, interpolated[i, ws_indxs])
    # ax.plot(DG.sizes[i] * 1e4, em_3, 'o', color="black", label="inter")

    # em_5 = np.interp(args.wave, waves, DG.emission_5[i, ws_indxs])
    # ax.plot(DG.sizes[i] * 1e4, em_5, 'o', color="green", label="5")

    # ax.set_xlabel(r"$a$ [$\mu m$]")
    # if args.obsdata != "none":
    # ax.set_ylabel(f"S({args.wave})/A(V)")
    # else:
    # ax.set_ylabel(f"S({args.wave})/N(HI)")
    # ax.set_xscale("log")
    # ax.set_yscale("log")

    size = DG.sizes[args.grainsize] * 1e4
    em_05 = np.interp(args.wave, waves, DG.emission_05[args.grainsize, ws_indxs])
    ax.plot(0.5, em_05, "o", color="blue")

    em_1 = np.interp(args.wave, waves, DG.emission_1[args.grainsize, ws_indxs])
    ax.plot(1, em_1, "o", color="blue")

    em_2 = np.interp(args.wave, waves, DG.emission_2[args.grainsize, ws_indxs])
    ax.plot(2, em_2, "o", color="blue")

    em_5 = np.interp(args.wave, waves, DG.emission_5[args.grainsize, ws_indxs])
    ax.plot(5, em_5, "o", color="blue")

    em_10 = np.interp(args.wave, waves, DG.emission_10[args.grainsize, ws_indxs])
    ax.plot(10, em_10, "o", color="blue")

    ax.set_title(f"grainsize: {size} $\u03bcm$, wavelength: {args.wave} $\u03bcm$")
    ax.set_xlabel("ISRF strength")
    ax.set_ylabel("Emission")

    # tot_em_05 = np.sum(DG.emission_05, axis=0)
    # tot_em_1 = np.sum(DG.emission_1, axis=0)
    # tot_em_2 = np.sum(DG.emission_2, axis=0)
    # tot_em_5 = np.sum(DG.emission_5, axis=0)
    # tot_em_10 = np.sum(DG.emission_10, axis=0)
    # interpolated_array = np.sum(interpolated, axis=0)

    # ax.plot(DG.wavelengths_emission[ws_indxs], DG.emission_05[args.wave, ws_indxs], label="0.5")
    # ax.plot(waves, tot_em_1[ws_indxs], label="1")
    # ax.plot(waves, interpolated_array[ws_indxs], label="inter")
    # ax.plot(waves, tot_em_2[ws_indxs], label="2")
    # ax.plot(waves, tot_em_5[ws_indxs], label="5")
    # ax.plot(DG.wavelengths_emission[ews_indxs], DG.emission_10[args.wave, ews_indxs], label="10")

    # ax.set_xlim(1, 10**4)
    # ax.set_ylim(10**(-20), 1)
    # ax.set_xscale("log")
    # ax.set_yscale("log")

    plt.tight_layout()

    # show or save
    basename = "DustGrains_diag_%s" % (args.composition)
    if args.png:
        fig.savefig(basename + ".png")
    elif args.eps:
        fig.savefig(basename + ".eps")
    elif args.pdf:
        fig.savefig(basename + ".pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()
