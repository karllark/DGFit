############
Plot Results
############

A number of programs exist to plot the inputs and results of the fitting.

Results
=======

To plot the results in filename based on the obsdata, use the command 

    plot_dgfit filename obsdata

For example, to plot the default run results from the optimizer (= max prob):

    plot_dgfit dgfit_test_WD_best_optimizer.fits obsdata

.. plot::

    from dgfit.plotting.plot_dgfit import (plot_dgfit_sizedist,
                                            plot_dgfit_extinction,
                                            plot_dgfit_emission,
                                            plot_dgfit_albedo,
                                            plot_dgfit_g,
                                            plot_dgfit_abundances)
    from dgfit.obsdata import ObsData
    import numpy as np
    import matplotlib.pyplot as pyplot
    import matplotlib
    from astropy.io import fits

    OD = ObsData("mw_rv31_obs.dat", path = "../../dgfit/data/mw_rv31/")
    hdulist = fits.open("../../dgfit/data/mw_rv31/dgfit_test_WD_sizedist_best_optimizer.fits")

    fontsize = 16
    font = {"size": fontsize}
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = pyplot.subplots(ncols=3, nrows=2, figsize=(15, 10))

    plot_dgfit_sizedist(ax[0, 0], hdulist, fontsize=fontsize, plegend=True)
    plot_dgfit_abundances(
        ax[0, 1],
        hdulist["ABUNDANCES"],
        OD,
        fontsize=fontsize,
        color="r",
        plegend=True,
        plabel="Final",
    )
    plot_dgfit_extinction(ax[1, 0], hdulist["EXTINCTION"], OD, fontsize=fontsize)
    plot_dgfit_emission(ax[0, 2], hdulist["EMISSION"], OD, fontsize=fontsize)
    plot_dgfit_albedo(ax[1, 1], hdulist["ALBEDO"], OD, fontsize=fontsize)
    plot_dgfit_g(ax[1, 2], hdulist["G"], OD, fontsize=fontsize)

    ax[0, 0].set_ylim(1e-40, 3e-27)
    pyplot.tight_layout()

The figure should look like this.

To include the starting points

    plot_dgfit dgfit_test_WD_best_optimizer.fits obsdata --start

To plot the data of the used dustgrains (default is astro-silicates), use command

    plot_dustgrains

.. plot::

    from dgfit.dustgrains import DustGrains
    import numpy as np
    import matplotlib.pyplot as pyplot
    import matplotlib
    import colorsys

    DG = DustGrains()
    DG.from_files("astro-silicates", "../../dgfit/data/indiv_grain/")

    fontsize = 12
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))

    ws_indxs = np.argsort(DG.wavelengths)
    ews_indxs = np.argsort(DG.wavelengths_emission)
    waves = DG.wavelengths[ws_indxs]
    colors = []
    labels = []
    for i in range(DG.n_sizes):
        pcolor = colorsys.hsv_to_rgb(float(i) / DG.n_sizes / (1.1), 1, 1)
        colors.append(pcolor)
        labels.append(DG.sizes[i])

        ax[0, 0].plot(waves, DG.cabs[i, ws_indxs], color=pcolor)
        ax[0, 0].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[0, 0].set_ylabel("C(abs)")
        ax[0, 0].set_xscale("log")
        ax[0, 0].set_yscale("log")

        ax[0, 1].plot(waves, DG.csca[i, ws_indxs], color=pcolor)
        ax[0, 1].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[0, 1].set_ylabel("C(sca)")
        ax[0, 1].set_xscale("log")
        ax[0, 1].set_yscale("log")

        ax[0, 2].plot(waves, DG.cext[i, ws_indxs], color=pcolor)
        ax[0, 2].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[0, 2].set_ylabel("C(ext)")
        ax[0, 2].set_xscale("log")
        ax[0, 2].set_yscale("log")

        ax[1, 0].plot(
            DG.wavelengths_scat_a,
            DG.scat_a_csca[i, :] / DG.scat_a_cext[i, :],
            "o",
            color=pcolor,
        )
        ax[1, 0].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[1, 0].set_ylabel("albedo")
        ax[1, 0].set_xscale("log")

        ax[1, 1].plot(DG.wavelengths_scat_g, DG.scat_g[i, :], "o", color=pcolor)
        ax[1, 1].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[1, 1].set_ylabel("g")
        ax[1, 1].set_xscale("log")

        ax[1, 2].plot(
            DG.wavelengths_emission[ews_indxs], DG.emission[i, ews_indxs], color=pcolor
        )
        ax[1, 2].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[1, 2].set_ylabel("Emission")
        ax[1, 2].set_xscale("log")
        ax[1, 2].set_yscale("log")
        cur_ylim = ax[1, 2].get_ylim()
        ax[1, 2].set_ylim([1e-23, 1e-0])

    ax[0, 1].set_title("astro-silicates")
    fig.legend(labels, title="Grainsizes [$m$]", loc="lower center", bbox_to_anchor=(0.5, 0), ncol=DG.n_sizes/3)


    plt.tight_layout(rect=[0, 0.14, 1, 1])

The figure will look like this.

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

    plot_dustgrains -c=<possible>

To transform the particles to the observed data grids:

    plot_dustgrains --obsdata obsdata

To see the options for saving the plots, use

    plot_dustgrains -help

To see an overview of the observed data used, use

    plot_obsdata filename

.. plot::

    import numpy as np
    import matplotlib.pyplot as pyplot
    import matplotlib
    from astropy.table import Table

    from dgfit.obsdata import ObsData

    OD = ObsData("mw_rv31_obs.dat", path = "../../dgfit/data/mw_rv31/")

    fontsize = 16
    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=2)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("ytick.major", width=2)

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))

    ax[0, 0].errorbar(
        OD.ext_waves, OD.ext_alnhi, yerr=OD.ext_alnhi_unc, fmt="o", label="Extinction"
    )
    ax[0, 0].set_xlabel(r"$\lambda [\mu m]$")
    ax[0, 0].set_ylabel(r"$A(\lambda)/N(HI)$")
    ax[0, 0].set_xscale("log")
    ax[0, 0].set_xlim(0.085, 3.0)
    ax[0, 0].legend()

    n_atoms = len(OD.abundance)
    aindxs = np.arange(n_atoms)
    width = 0.5
    atomnames = sorted(list(OD.abundance.keys()))

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

    ax[1, 0].set_ylabel(r"$N(X)/[10^6 N(HI)]$", fontsize=fontsize)
    ax[1, 0].set_xticks(aindxs + (0.75 * width))
    ax[1, 0].set_xticklabels(atomnames)
    ax[1, 0].legend(loc=2)

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

This will look like this, depending on what obsdata is used.

To add the ISRF plot (if available)

    plot_obsdata filename --ISRF ISRFdatafile

This ISRF plot will pop up in the middle plot of the lower row.
