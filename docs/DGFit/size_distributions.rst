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

WD
==

`Weingartner & Draine (2001)
<https://ui.adsabs.harvard.edu//#abs/2001ApJ...548..296W/abstract>`_
size distributions ...
