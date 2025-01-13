from __future__ import print_function

import argparse
import importlib.resources as importlib_resources

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy.table import Table

from dgfit.obsdata import ObsData


def main():

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "obsdata", type=str, default="none", help="give the file with the observed data"
    )
    parser.add_argument(
        "--ISRF", type=str, default="none", help="Add the ISFR file to plot"
    )
    parser.add_argument(
        "--units",
        default="AV",
        choices=["AV", "NHI"],
        help="Choose in what units the plots are",
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    OD = ObsData(args.obsdata)

    plot(OD, args.ISRF, args.units, args.png, args.eps, args.pdf)


def plot(OD, ISRF="none", units="AV", png=False, eps=False, pdf=False):

    # setup the plots
    fontsize = 16
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))

    n_atoms = len(OD.abundance)
    aindxs = np.arange(n_atoms)
    width = 0.5
    atomnames = sorted(list(OD.abundance.keys()))

    if units == "AV":
        ax[0, 0].errorbar(
            OD.ext_waves, OD.ext_alav, yerr=OD.ext_alav_unc, fmt="o", label="Extinction"
        )
        ax[0, 0].set_ylabel(r"$A(\lambda)/A(V)$")

        ax[1, 0].bar(
            aindxs + 0.25 * width,
            [OD.total_abundance_av[x][0] for x in atomnames],
            width,
            color="g",
            alpha=0.25,
            label="gas+dust",
        )
        ax[1, 0].errorbar(
            aindxs + 0.75 * width,
            [OD.abundance_av[x][0] for x in atomnames],
            yerr=[OD.abundance_av[x][1] for x in atomnames],
            fmt="o",
            label="dust",
        )
        ax[1, 0].set_ylabel(r"$N(X)/A(V)$", fontsize=fontsize)

        if OD.fit_ir_emission:
            ax[0, 1].errorbar(
                OD.ir_emission_waves,
                OD.ir_emission_av,
                yerr=OD.ir_emission_av_unc,
                fmt="o",
                label="Emission",
            )
            ax[0, 1].set_xlabel(r"$\lambda [\mu m]$")
            ax[0, 1].set_ylabel(r"$S$ $[MJy$ $sr^{-1}$ $A(V)^{-1}]$")
            ax[0, 1].set_xscale("log")
            ax[0, 1].set_xlim(1.0, 1.5e4)
            ax[0, 1].set_yscale("log")
            ax[0, 1].legend(loc=2)

    elif units == "NHI":
        ax[0, 0].errorbar(
            OD.ext_waves,
            OD.ext_alnhi,
            yerr=OD.ext_alnhi_unc,
            fmt="o",
            label="Extinction",
        )
        ax[0, 0].set_ylabel(r"$A(\lambda)/N(HI)$")

        ax[1, 0].bar(
            aindxs + 0.25 * width,
            [OD.total_abundance[x][0] for x in atomnames],
            width,
            color="g",
            alpha=0.25,
            label="gas+dust",
        )
        ax[1, 0].errorbar(
            aindxs + 0.75 * width,
            [OD.abundance[x][0] for x in atomnames],
            yerr=[OD.abundance[x][1] for x in atomnames],
            fmt="o",
            label="dust",
        )
        ax[1, 0].set_ylabel(r"$N(X)/[10^6N(HI)]$", fontsize=fontsize)

        if OD.fit_ir_emission:
            ax[0, 1].errorbar(
                OD.ir_emission_waves,
                OD.ir_emission,
                yerr=OD.ir_emission_unc,
                fmt="o",
                label="Emission",
            )
            ax[0, 1].set_xlabel(r"$\lambda [\mu m]$")
            ax[0, 1].set_ylabel(r"$S$ $[MJy$ $sr^{-1}$ $N(HI)^{-1}]$")
            ax[0, 1].set_xscale("log")
            ax[0, 1].set_xlim(1.0, 1.5e4)
            ax[0, 1].set_yscale("log")
            ax[0, 1].legend(loc=2)

    ax[0, 0].set_xlabel(r"$\lambda [\mu m]$")
    ax[0, 0].set_xscale("log")
    ax[0, 0].set_xlim(0.085, 3.0)
    ax[0, 0].legend()

    ax[1, 0].set_xticks(aindxs + (0.75 * width))
    ax[1, 0].set_xticklabels(atomnames)
    ax[1, 0].legend(loc=2)

    if ISRF != "none":
        ref = importlib_resources.files("dgfit") / ISRF
        with importlib_resources.as_file(ref) as data_path:
            t = Table.read(str(data_path), format="ascii.commented_header")
        ax[1, 1].plot(t["wave"], t["ISRF"], "-", label="ISRF")
        ax[1, 1].set_xlabel(r"$\lambda [\mu m]$")
        ax[1, 1].set_ylabel(r"ISRF [$ergs$ $cm^{-3}$ $s^{-1}$ $sr^{-1}$]")
        ax[1, 1].set_xscale("log")
        ax[1, 1].set_yscale("log")
        ax[1, 1].set_xlim(0.09, 1e1)
        ax[1, 1].set_ylim(1e-2, 1e2)
        ax[1, 1].legend()

    if OD.fit_scat_a:
        ax[0, 2].errorbar(
            OD.scat_a_waves,
            OD.scat_albedo,
            yerr=OD.scat_albedo_unc,
            fmt="o",
            label="albedo",
        )
        ax[0, 2].set_xlabel(r"$\lambda [\mu m]$")
        ax[0, 2].set_ylabel(r"$a$")
        ax[0, 2].set_xscale("log")
        ax[0, 2].set_xlim(0.085, 3.0)
        ax[0, 2].set_ylim(0.0, 1.0)
        ax[0, 2].legend()

    if OD.fit_scat_g:
        ax[1, 2].errorbar(
            OD.scat_g_waves,
            OD.scat_g,
            yerr=OD.scat_g_unc,
            fmt="o",
            label=r"$g = < \mathrm{cos} (\theta) >$",
        )
        ax[1, 2].set_xlabel(r"$\lambda [\mu m]$")
        ax[1, 2].set_ylabel(r"$g$")
        ax[1, 2].set_xscale("log")
        ax[1, 2].set_xlim(0.085, 3.0)
        ax[1, 2].set_ylim(0.0, 1.0)
        ax[1, 2].legend()

    plt.tight_layout()

    # show or save
    basename = "ObsData_MW_Diffuse"
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
