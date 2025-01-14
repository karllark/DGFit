#####################
Dust grain properties
#####################

There are several plots that show the effective grain properties of the used dustgrains.
In all plots, it is possible to transform the particles and data to the observed data that was used.

To plot the extinction coëfficients, albedo, scattering phase function and emission of the used dustgrains (default is astro-silicates) in function of the grain size, use command

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

    from dgfit.dustmodel import DustModel
    import dgfit.plotting.plot_effsize
    import importlib.resources as importlib_resources

    compnames = ["astro-silicates", "astro-carbonaceous"]
    ref = importlib_resources.files("dgfit") / "data"
    with importlib_resources.as_file(ref) as data_path:
        dustmodel_full = DustModel(
            componentnames=compnames,
            path=str(data_path) + "/indiv_grain/",
            every_nth=1,
        )
    dustmodel = WDDustModel(dustmodel=dustmodel_full)

    dgfit.plotting.plot_effsize.plot(dustmodel)


