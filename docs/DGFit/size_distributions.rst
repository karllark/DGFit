##################
Size Distributions
##################

DGFit has different functional forms for size distributions that can be
used when fitting observed data.

Arbitrary
=========

This is the most flexible size distribution where each size bin is a
separate parameter.  The total number of parameters is equal to the
sum of the number of size distributions for all the different grain
compositions.

MRN
===

The commonly used `Mathis, Rumpl, & Nordsieck (1977)
<https://ui.adsabs.harvard.edu//#abs/1977ApJ...217..425M/abstract>`_
size distribution is a powerlaw of the form:

.. math::
  \begin{eqnarray}
    w(a) & = C a ^{-\alpha} \quad\quad & a_{min} \leq a \leq a_{max} \\
    w(a) & = 0 \quad\quad\quad & a < a_{min} \quad \mathrm{or} \quad a > a_{max} \\
  \end{eqnarray}

where :math:`C` is the amplitude, :math:`a` is the grain radius, and
:math:`\alpha` is the powerlaw index.

.. plot::
    import numpy as np
    import matplotlib.pyplot as plt

    from DGFit.DustModel import MRNDustModel

    fig, ax = plt.subplots()

    # simple intialization to get the MRN size distribution function
    dmodel = MRNDustModel()

    # size distributions in cm
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)*1e-8

    # default size distribution parameters
    cparams = [1e-25, 3.5, 1e-7, 1e-3]

    sizedist = dmodel.compute_size_dist(a, cparams)

    ax.plot(a, sizedist)

    # ax.set_xlabel('$x$ [$\mu m^{-1}$]')
    # ax.set_ylabel('$A(x)/A(V)$')

    # ax.legend(loc='best')
    plt.tight_layout()
    plt.show()

WD
==

`Weingartner & Draine (2001)
<https://ui.adsabs.harvard.edu//#abs/2001ApJ...548..296W/abstract>`_
size distributions are defined in equations 4 (carbonaceous grains)
& 5 (silicate grains) of the paper.
