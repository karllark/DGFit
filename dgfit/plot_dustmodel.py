import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from dgfit.obsdata import ObsData
from dgfit.dustmodel import DustModel

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--obsdata", help="transform to observed data grids", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    path = "dgfit/data/mw_rv31"
    OD = ObsData(
        [
            "%s/MW_diffuse_Gordon09_band_ext.dat" % path,
            "%s/MW_diffuse_Gordon09_iue_ext.dat" % path,
            "%s/MW_diffuse_Gordon09_fuse_ext.dat" % path,
        ],
        "%s/MW_diffuse_Gordon09_avnhi.dat" % path,
        "%s/MW_diffuse_Jenkins09_abundances.dat" % path,
        "%s/MW_diffuse_Compiegne11_ir_emission.dat" % path,
        "%s/dust_scat.dat" % path,
        ext_tags=["band", "iue", "fuse"],
        scat_path="%s/Scat_Data/" % path,
    )

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    # dustmodel = DustModel(['astro-silicates','astro-graphite'])
    DM = DustModel()
    DM.predict_full_grid(
        ["astro-silicates", "astro-carbonaceous"], path="dgfit/data/indiv_grain/"
    )

    if args.obsdata:
        DM_obs = DustModel()
        DM_obs.predict_observed_data(DM, OD)
        DM = DM_obs

    results = DM.eff_grain_props(OD)
    cabs = results["cabs"]
    csca = results["csca"]
    natoms = results["natoms"]
    emission = results["emission"]
    albedo = results["albedo"]
    g = results["g"]

    fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(14, 10))

    # plot the total results
    ax[1, 0].plot(DM.components[0].wavelengths, cabs + csca, "k-")
    ax[1, 0].set_xscale("log")
    ax[1, 0].set_yscale("log")
    ax[1, 0].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[1, 0].set_ylabel(r"C(ext)")
    # ax[1,0].set_xlim(1e-2,1e0)
    # ax[1,0].set_ylim(1e3,1e5)

    ax[1, 1].plot(DM.components[0].wavelengths_emission, emission, "k-")
    ax[1, 1].set_xscale("log")
    ax[1, 1].set_yscale("log")
    ax[1, 1].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[1, 1].set_ylabel(r"S")
    # ax[1,1].set_xlim(1e0,1e4)
    (gindxs,) = np.where(DM.components[0].wavelengths > 1e0)
    # ax[1,1].set_ylim(min(emission[gindxs]), max(emission[gindxs]))

    ax[1, 2].plot(DM.components[0].wavelengths_scat_a, albedo, "k-")
    ax[1, 2].set_xscale("log")
    ax[1, 2].set_yscale("linear")
    ax[1, 2].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[1, 2].set_ylabel(r"$a$")

    ax[2, 2].plot(DM.components[0].wavelengths_scat_g, g, "k-")
    ax[2, 2].set_xscale("log")
    ax[2, 2].set_yscale("linear")
    ax[2, 2].set_xlabel(r"$\lambda$ [$\mu m$]")
    ax[2, 2].set_ylabel(r"$g$")

    # plot the size distributions and component results
    for component in DM.components:
        ax[0, 0].plot(component.sizes, component.size_dist, "-", label=component.name)
        ax[0, 1].plot(
            component.sizes,
            np.power(component.sizes, 4.0) * component.size_dist,
            "-",
            label=component.name,
        )
        ax[0, 2].plot(
            component.sizes,
            np.power(component.sizes, 3.0) * component.size_dist,
            "-",
            label=component.name,
        )

        cresults = component.eff_grain_props(OD)
        ax[1, 0].plot(component.wavelengths, cresults["cabs"] + cresults["csca"])
        ax[1, 1].plot(component.wavelengths_emission, cresults["emission"])
        ax[1, 2].plot(component.wavelengths_scat_a, cresults["albedo"])
        ax[2, 2].plot(component.wavelengths_scat_g, cresults["g"])

    ax[0, 0].set_xscale("log")
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlabel(r"$a$ [$cm$]")
    ax[0, 0].set_ylabel(r"$N(a)$")

    ax[0, 1].set_xscale("log")
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_xlabel(r"$a$ [$cm$]")
    ax[0, 1].set_ylabel(r"$a^4 N(a)$")

    ax[0, 2].set_xscale("log")
    ax[0, 2].set_yscale("log")
    ax[0, 2].set_xlabel(r"$a$ [$cm$]")
    ax[0, 2].set_ylabel(r"$a^3 N(a)$")

    ax[0, 0].legend()

    plt.tight_layout()

    # show or save
    basename = "ObsData_MW_Diffuse"
    if args.png:
        fig.savefig(basename + ".png")
    elif args.eps:
        fig.savefig(basename + ".eps")
    elif args.pdf:
        fig.savefig(basename + ".pdf")
    else:
        plt.show()
