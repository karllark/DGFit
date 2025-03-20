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
        "dustproperty",
        help="What dust property needs to be shown for both model and data",
        choices=["emission", "extinction", "albedo", "g"],
    )
    parser.add_argument(
        "--start", help="include the starting model", action="store_true"
    )
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
    colors = ["r", "b", "g", "c"]
    ltype = "-"

    # open the DGFit results
    hdulist = fits.open(args.filename)

    comps = []
    for i in range(hdulist[0].header["NCOMPS"]):
        hdu = hdulist[i + 1]
        comps.append(hdu.header["EXTNAME"])

    hdu = hdulist[args.dustproperty.upper()]

    # get the observed data
    OD = ObsData(args.obsfile)

    if args.dustproperty == "emission":
        waves = OD.ir_emission_waves
        data = OD.ir_emission_av
        data_unc = OD.ir_emission_av_unc
        data_name = "EMIS"
        ylabel = r"$S$ $[MJy$ $sr^{-1}$ $A(V)^{-1}]$"
        logscale = True
        ylim = False

    elif args.dustproperty == "extinction":
        waves = OD.ext_waves
        data = OD.ext_alav
        data_unc = OD.ext_alav_unc
        data_name = "EXT"
        ylabel = r"$A(\lambda)/A(V)$"
        logscale = True
        ylim = False

    elif args.dustproperty == "albedo":
        waves = OD.scat_a_waves
        data = OD.scat_albedo
        data_unc = OD.scat_albedo_unc
        data_name = args.dustproperty.upper()
        ylabel = "Albedo"
        logscale = False
        ylim = True

    elif args.dustproperty == "g":
        waves = OD.scat_g_waves
        data = OD.scat_g
        data_unc = OD.scat_g_unc
        data_name = args.dustproperty.upper()
        ylabel = "g"
        logscale = False
        ylim = True

    ax1.plot(hdu.data["WAVE"], hdu.data[data_name], colors[0] + ltype, label="Total")
    yrange = get_krange(hdu.data[data_name])
    for i in range(len(hdu.data.names) - 2):
        ax1.plot(
            hdu.data["WAVE"],
            hdu.data[data_name + str(i + 1)],
            colors[i + 1] + ltype,
            label=comps[i],
        )
        yrange = get_krange(hdu.data[data_name + str(i + 1)], in_range=yrange)

    ax1.errorbar(
        waves,
        data,
        data_unc,
        fmt="ko",
        label="Observed",
        capsize=3,
    )

    if logscale:
        ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$\lambda [\mu m]$", fontsize=fontsize)
    ax1.set_ylabel(ylabel, fontsize=fontsize)
    ax1.legend()
    ax1.set_xlim(get_krange(hdu.data["WAVE"], logaxis=True))
    if ylim:
        ax1.set_ylim([0.0, 1.0])

    residuals = (hdu.data[data_name] - data) / data
    residuals *= 100
    unc = data_unc / data
    unc *= 100
    ax2.errorbar(
        hdu.data["WAVE"],
        residuals,
        yerr=unc,
        fmt="o",
        color="black",
        capsize=3,
    )
    ax2.axhline(0, color="red", linestyle="--", linewidth=1)
    ax2.set_xlabel(r"$\lambda [\mu m]$", fontsize=fontsize)
    ax2.set_ylabel("Residuals (%)", fontsize=fontsize)

    pyplot.tight_layout()

    # show or save
    basename = args.filename
    basename.replace(".fits", "")
    basename += "_" + args.dustproperty
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
