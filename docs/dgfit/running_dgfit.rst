################
How to run DGFit
################

The run DGFit, the dgfit package is assumed to be installed.  See <dgfit/install.rst>.

Commandline
===========

The fitting against the observed data file is run by using

.. code-block:: console

    run_dgfit obsdata

To run the default MW Rv=3.1 case quickly to test the installation, use

.. code-block:: console

    run_dgfit obdata --fast

The default run is done assuming the WD size distributions and the base filename
of `dgfit_test_WD`.  This default run will be done with a coarse set of dust
grain sizes (every 2nd possible).

To run the with all the dust grains sizes possible

.. code-block:: console

    run_dgfit obsdata --everynth=1

To run with different size distributions (<possible> = MRN, WD, bins), use:

.. code-block:: console

    run_dgfit obsdata --sizedisttype=<possible>

For all the `run_dgfit` options use:

.. code-block:: console
    
    run_dgfit --help