import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from astropy.io import fits

from dgfit.obsdata import ObsData


def get_krange(x, logaxis=False, in_range=[0]):
    prange = np.array([0.0, 0.0])
    if logaxis:
        gindxs = x > 0
        min_x = np.amin(x[gindxs])
        max_x = np.amax(x[gindxs])
    else:
        min_x = np.amin(x)
        max_x = np.amax(x)

    if logaxis:
        max_x = np.log10(max_x)
        if min_x <= 0.0:
            min_x = max_x - 10.0
        else:
            min_x = np.log10(min_x)

    delta = max_x - min_x
    prange[0] = min_x - 0.1 * delta
    prange[1] = max_x + 0.1 * delta

    if logaxis:
        prange = np.power(10.0, prange)

    if len(in_range) > 1:
        prange[0] = np.minimum(prange[0], in_range[0])
        prange[1] = np.maximum(prange[1], in_range[1])

    return prange


def main():
    # commandline parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "filename",
        help=(
            "file with the dust model details "
            + "(size distribution, extinction, etc.)"
        ),
    )
    parser.add_argument(
        "obsfile", help="Data file giving the observational data that was fit"
    )
    parser.add_argument(
        "--start", help="include the starting model", action="store_true"
    )
    parser.add_argument("--smc", help="use an SMC sightline", action="store_true")
    parser.add_argument(
        "-p", "--png", help="save figure as a png file", action="store_true"
    )
    parser.add_argument(
        "-e", "--eps", help="save figure as an eps file", action="store_true"
    )
    parser.add_argument("-pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # setup the plots
    fontsize = 16
    font = {"size": fontsize}

    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(
        ncols=1,
        nrows=2,
        figsize=(15, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    ax1 = ax[0]
    ax2 = ax[1]
    colors = ["r", "b", "g"]
    ltype = "-"

    # open the DGFit results
    hdulist = fits.open(args.filename)

    comps = []
    for i in range(hdulist[0].header["NCOMPS"]):
        hdu = hdulist[i + 1]
        comps.append(hdu.header["EXTNAME"])

    hdu = hdulist["EMISSION"]

    # get the observed data
    OD = ObsData(args.obsfile)

    ax1.plot(hdu.data["WAVE"], hdu.data["EMIS"], colors[0] + ltype, label="Total")
    yrange = get_krange(hdu.data["EMIS"], logaxis=True)
    for i in range(len(hdu.data.names) - 2):
        ax1.plot(
            hdu.data["WAVE"],
            hdu.data["EMIS" + str(i + 1)],
            colors[i + 1] + ltype,
            label=comps[i],
        )
        yrange = get_krange(
            hdu.data["EMIS" + str(i + 1)], logaxis=True, in_range=yrange
        )

    ax1.errorbar(
        OD.ir_emission_waves,
        OD.ir_emission_av,
        yerr=OD.ir_emission_av_unc,
        fmt="ko",
        label="Observed",
        capsize=4,
    )
    yrange = get_krange(OD.ir_emission_av, logaxis=True, in_range=yrange)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend()
    ax1.set_xlabel(r"$\lambda [\mu m]$", fontsize=fontsize)
    ax1.set_ylabel(r"$S$ $[MJy$ $sr^{-1}$ $A(V)^{-1}]$", fontsize=fontsize)
    ax1.set_xlim(get_krange(hdu.data["WAVE"], logaxis=True))
    ax1.set_ylim(yrange)

    residuals = (hdu.data["EMIS"] - OD.ir_emission_av) / OD.ir_emission_av
    unc = OD.ir_emission_av_unc / OD.ir_emission_av
    ax2.errorbar(hdu.data["WAVE"], residuals, yerr=unc, fmt="o", color='black', capsize=4)
    ax2.axhline(0, color="red", linestyle="--", linewidth=1)
    ax2.set_xlabel(r"$\lambda [\mu m]$", fontsize=fontsize)
    ax2.set_ylabel("Residuals (%)", fontsize=fontsize)

    pyplot.tight_layout()

    # show or save
    basename = args.filename
    basename.replace(".fits", "")
    if args.png:
        fig.savefig(basename + ".png")
    elif args.eps:
        fig.savefig(basename + ".eps")
    elif args.pdf:
        fig.savefig(basename + ".pdf")
    else:
        pyplot.show()


if __name__ == "__main__":
    main()
