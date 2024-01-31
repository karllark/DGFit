############
Plot Results
############

A number of programs exist to plot the inputs and results of the fitting.

Results
=======

To plot the results, use the command 

    plot_dgfit filename

For example, to plot the default run results from the optimizer (= max prob):

    plot_dgfit dgfit_test_WD_best_optimizer.fits

To include the starting points

    plot_dgfit dgfit_test_WD_best_optimizer.fits --start