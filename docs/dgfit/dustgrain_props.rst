#####################
Dust grain properties
#####################

There are several plots that show the effective grain properties of the used dustgrains.
In all plots, it is possible to transform the particles and data to the observed data that was used.

General properties
------------------

To plot the extinction coefficients, albedo, scattering phase function and emission of the used dustgrains (default is astro-silicates) in function of the grain size, use command

.. code-block:: console

    dgplot_dustgrains

.. plot::

    from dgfit.dustgrains import DustGrains
    from dgfit.plotting.plot_dustgrains import plot

    DG = DustGrains()
    DG.from_files("astro-silicates", "../../dgfit/data/indiv_grain/")

    plot(DG, 'astro-silicates', 1)

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

.. code-block:: console

    dgplot_dustgrains -c <possible>

To transform the particles to the observed data grids and see the data of the dustgrains for the observed dustmodel, use

.. code-block:: console

    dgplot_dustgrains --obsdata obsdata

To see the options for saving the plots (for this and future plots), use

.. code-block:: console

    dgplot_dustgrains --help

Effective size
--------------

To plot the average dustgrain size in function of wavelength for both extinction and emission, use

.. code-block:: console

    dgplot_effsize

.. plot::

    from dgfit.dustmodel import DustModel, WDDustModel
    from dgfit.dustgrains import DustGrains
    from dgfit.plotting.plot_effsize import plot
    
    compnames = ["astro-silicates", "astro-carbonaceous"]

    DustM = DustModel(componentnames=compnames, path="../../dgfit/data/indiv_grain/", dustmodel=None, obsdata=None, every_nth=2)
    model = WDDustModel(dustmodel=DustM)

    plot(model)

Single grain properties
-----------------------

Another plot that is available is one that shows the single grain properties. It is devided into four subplots.
The first one shows the extinction for a chosen wavelength and composition in function of grain size.
The second one shows the same thing but for the emission.
The third and fourth subplot show the relative absorpiton and scattering coefficient for a chosen wavelength in function of grain size.
To make this plot, use

.. code-block:: console

    dgplot_singlegrain_props

.. plot::

    from dgfit.dustgrains import DustGrains
    from dgfit.plotting.plot_singlegrain_props import plot

    DG = DustGrains()
    DG.from_files("astro-silicates", "../../dgfit/data/indiv_grain/", every_nth=1)

    plot(DG, 10, 'astro-silicates', 1)

The default wavelength for these plots is 0.1 micron. 
To choose another wavelength (with the wavelength given in micron), use

.. code-block:: console

    dgplot_singlegrain_props --wave wavelength

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

.. code-block:: console

    dgplot_singlegrain_props -c <possible>

To transform the particles to the observed data grids and see the data of the dustgrains for the observed dustmodel, use

.. code-block:: console

    dgplot_singlegrain_props --obsdata obsdata

R(V)
----

The last available plotting options shows the extinction-R(V) relation for all different grain sizes and a chosen wavelength and composition.
To make this plot, use

.. code-block:: console

    dgplot_rv

.. plot::

    from dgfit.dustgrains import DustGrains
    from dgfit.plotting.plot_rv import plot

    DG = DustGrains()
    DG.from_files("astro-silicates", "../../dgfit/data/indiv_grain/", every_nth=1)

    plot(DG, 10, 'astro-silicates')

The default wavelength for these plots is 0.1 micron. 
To choose another wavelength (with the wavelength given in micron), use

.. code-block:: console

    dgplot_rv --wave wavelength

To see other dustgrains (<possible> = astro-silicates, astro-carbonaceous, astro-graphite, astro-PAH-ionized and astro-PAH-neutral), use

.. code-block:: console

    dgplot_rv -c <possible>

To transform the particles to the observed data grids and see the data of the dustgrains for the observed dustmodel, use

.. code-block:: console
    
    dgplot_rv --obsdata obsdata