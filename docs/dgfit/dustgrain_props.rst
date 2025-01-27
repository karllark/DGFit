#####################
Dust grain properties
#####################

There are several plots that show the effective grain properties of the used dustgrains.
In all plots, it is possible to transform the particles and data to the observed data that was used.

To plot the extinction coefficients, albedo, scattering phase function and emission of the used dustgrains (default is astro-silicates) in function of the grain size, use command

    dgplot_dustgrains

.. plot::

    from dgfit.dustgrains import DustGrains
    import dgfit.plotting.plot_dustgrains

    DG = DustGrains()
    DG.from_files("astro-silicates", "../../dgfit/data/indiv_grain/")

    dgfit.plotting.plot_dustgrains.plot(DG, 'astro-silicates')

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

    dgplot_dustgrains -c <possible>

To transform the particles to the observed data grids and see the data of the dustgrains for the observed dustmodel:

    dgplot_dustgrains --obsdata obsdata

To see the options for saving the plots (for this and future plots), use

    dgplot_dustgrains --help

To plot the average dustgrain size in function of wavelength for both extinction and emission, use

    dgplot_effsize

.. plot::

    from dgfit.dustmodel import DustModel, WDDustModel
    from dgfit.dustgrains import DustGrains
    import dgfit.plotting.plot_effsize
    
    compnames = ["astro-silicates", "astro-carbonaceous"]

    DM = DustModel(componentnames="astro-carbonaceous", path="../../dgfit/data/indiv_grain/", dustmodel=None, obsdata=None, every_nth=2)
    model = WDDustModel(DM)

    dgfit.plotting.plot_effsize.plot(model)



