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


.. .. figure:: plot_dgfit.png
..     :width: 50%
..     :align: center

The figure should look like this.

To include the starting points

    plot_dgfit dgfit_test_WD_best_optimizer.fits obsdata --start

To plot the data of the used dustgrains (default is astro-silicates), use command

    plot_dustgrains

.. figure:: plot_dustgrains.png
    :width: 50%
    :align: center

The figure will look like this.

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

    plot_dustgrains -c=<possible>

To transform the particles to the observed data grids:

    plot_dustgrains --obsdata obsdata

To see the options for saving the plots, use

    plot_dustgrains -help

To see an overview of the observed data used, use

    plot_obsdata filename

To add the ISRF plot (if available)

    plot_obsdata filename --ISRF ISRFdatafile

