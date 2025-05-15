import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from astropy.io import fits


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


def plot(
    ax,
    hdulist,
    colors=["b", "g", "c", "r"],
    markers=["x", "o", "s", "*"],
    fontsize=12,
    multa4=False,
    mass=False,
    plegend=True,
    ltype="-",
    alpha=1.0,
):
    if "DISTPUNC" in hdulist[1].data.names:
        plot_uncs = True
    else:
        plot_uncs = False

    yrange = [0]
    all_yvals = []
    for i in range(hdulist[0].header["NCOMPS"]):
        hdu = hdulist[i + 1]

        xvals = hdu.data["SIZE"] * 1e4
        yvals = hdu.data["DIST"]

        if np.sum(yvals) == 0:
            print(f"Composition {i} is zero")
            continue

        if plot_uncs:
            yvals_punc = hdu.data["DISTPUNC"]
            yvals_munc = hdu.data["DISTMUNC"]

        if multa4:
            xvals4 = hdu.data["SIZE"] ** 4
            yvals *= xvals4
            if plot_uncs:
                yvals_punc *= xvals4
                yvals_munc *= xvals4

        elif mass:
            xvals3 = hdu.data["SIZE"] ** 3
            yvals *= xvals3
            if plot_uncs:
                yvals_punc *= xvals3
                yvals_munc *= xvals3

        yrange = get_krange(yvals, logaxis=True, in_range=yrange)
        if plot_uncs:
            yrange = get_krange(yvals - yvals_munc, logaxis=True, in_range=yrange)
            yrange = get_krange(yvals + yvals_punc, logaxis=True, in_range=yrange)

        gindxs = yvals > 0
        all_yvals.append(np.max(yvals))

        ax.plot(
            xvals[gindxs],
            yvals[gindxs],
            colors[i] + ltype,
            marker=markers[i],
            markevery=3,
            label=hdu.header["EXTNAME"],
            alpha=alpha,
        )
        if plot_uncs:
            ax.errorbar(
                xvals[gindxs],
                yvals[gindxs],
                fmt=colors[i] + "o",
                yerr=[yvals_munc[gindxs], yvals_punc[gindxs]],
                alpha=alpha,
            )

    if multa4:
        ylabel = r"$a^4 N_d(a)/A(V)$"
    elif mass:
        ylabel = r"m(a)/A(V)"
    else:
        ylabel = r"$N_d(a)/A(V)$"

    ymax = max(all_yvals) * 100
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(get_krange(xvals, logaxis=True))
    if multa4:
        ax.set_ylim(1e-15, ymax)
    elif mass:
        ax.set_ylim(1e-5, ymax)
    else:
        ax.set_ylim(1e0, ymax)
    ax.set_xlabel(r"a $[\mu m]$", fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    if plegend:
        ax.legend()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filename",
        help=(
            "file with the dust model details "
            + "(size distribution, extinction, etc.)"
        ),
    )
    parser.add_argument(
        "--start", help="include the starting model", action="store_true"
    )
    parser.add_argument(
        "--startonly", help="include the starting model", action="store_true"
    )
    parser.add_argument(
        "--multa4", action="store_true", help="multiply the size distribution by a**4"
    )
    parser.add_argument(
        "--mass", action="store_true", help="Show the mass distribution"
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

    fontsize = 16
    font = {"size": fontsize}

    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(figsize=(15, 10))

    hdulist = fits.open(args.filename)
    if not args.startonly:
        plot(
            ax,
            hdulist,
            fontsize=fontsize,
            multa4=args.multa4,
            mass=args.mass,
            plegend=True,
        )

    if args.start or args.startonly:
        legend = False
        if args.startonly:
            legend = True
        if "best_fin" in args.filename:
            repstr = "best_fin"
        elif "fin" in args.filename:
            repstr = "fin"
        else:
            repstr = "best_optimizer"
        hdulist2 = fits.open(args.filename.replace(repstr, "start"))
        plot(
            ax,
            hdulist2,
            fontsize=fontsize,
            multa4=args.multa4,
            mass=args.mass,
            plegend=legend,
            ltype="--",
            alpha=0.50,
        )

    pyplot.tight_layout()

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
