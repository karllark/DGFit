################
How to run DGFit
################

The run DGFit, the dgfit package is assumed to be installed.  See <dgfit/install.rst>.

Commandline
===========

The fitting against the obsdata is run by using

    run_dgfit obsdata

To run the default MW Rv=3.1 case quickly to test the installation, use

    run_dgfit obdata --fast

The default run is done assuming the MRN size distributions and the base filename
of `dgfit_test_MRN`.  This default run will be done with a coarse set of dust
grain sizes (every 5th possible).

To run the with all the dust grains sizes possible

    run_dgfit obsdata --everynth=1

To run with different size distributions (<possible> = MRN, WD, bins), use:

    run_dgfit obsdata --sizedisttype=<possible>

For all the `run_dgfit` options use:

    run_dgfit --help