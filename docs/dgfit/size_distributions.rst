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
<https://ui.adsabs.harvard.edu/abs/1977ApJ...217..425M/abstract>`_
size distribution is a powerlaw of the form:

.. math::
  \frac{dn}{da} = \left\{
    \begin{array}{ll}
    C a ^{-\alpha} & a_{min} \leq a \leq a_{max} \\
    0 & a < a_{min} \\
    0 & a > a_{max} \\
    \end{array}
    \right.

where :math:`C` is the amplitude, :math:`a` is the grain radius,
:math:`\alpha` is the powerlaw index, and :math:`a_{min}`/:math:`a_{max}`
give the min/max grain radii.  A plot of an MRN size distribution
with :math:`C = 1\times 10^{-25}`, :math:`\alpha = 3.5`,
:math:`a_{min} = 3.5\times 10^{-4}~\mu m`,
and :math:`a_{max} = 2~\mu m` is shown below.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from dgfit.dustmodel import MRNDustModel

    fig, ax = plt.subplots()

    # simple intialization to get the MRN size distribution function
    dmodel = MRNDustModel()

    # size distributions in cm
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)*1e-8

    # default size distribution parameters
    cparams = [1e-25, 3.5, 1e-7, 1e-3]

    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, sizedist[indxs])

    ax.set_xlabel('$a$ [$\mu m$]')
    ax.set_ylabel(r'$n_H^{-1} \; dn/da$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.tight_layout()
    plt.show()

The same example size distributions are shown in a 2nd plot below
where the detailed wiggles in the overall powerlaw are emphasized
by multiplying by :math:`a^4`.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from dgfit.dustmodel import MRNDustModel

    fig, ax = plt.subplots()

    # simple intialization to get the MRN size distribution function
    dmodel = MRNDustModel()

    # size distributions in cm
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)*1e-8

    # default size distribution parameters
    cparams = [1e-25, 3.5, 1e-7, 1e-3]

    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, 1e29*(a[indxs]**4)*sizedist[indxs])

    ax.set_xlabel('$a$ [$\mu m$]')
    ax.set_ylabel(r'$10^{29} \; n_H^{-1} \; a^4 \; dn/da$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.tight_layout()
    plt.show()

WD
==

`Weingartner & Draine (2001)
<https://ui.adsabs.harvard.edu/abs/2001ApJ...548..296W/abstract>`_
size distributions are the referenced paper.  Combining the equations from
the paper results in a general form that applies to both carbonaceous and
silicate grains.  In the equation below, the first term provides a
curved power law, the 2nd term :math:`D(a)` is the
sum of two log-normal functions (only used for carbonaceous grains), and all is
multiplied by :math:`G(a)` that results in an exponential cutoff at
large grain radii.

.. math::
    \frac{dn}{da} = \left[ \frac{C}{a} \left( \frac{a}{a_t} \right)^\alpha
       F(a; \beta, \alpha) + D(a) \right] G(a)

where

.. math::
   F(a; \beta, \alpha) = \left\{
     \begin{array}{ll}
     1 + \beta a /a_t, & \beta \geq 0 \\
     (1 - \beta a/a_t)^{-1} & \beta < 0
     \end{array}
     \right.

and

.. math::
   G(a) = \left\{
     \begin{array}{ll}
     1, & 0.35~\mathrm{nm} < a < a_t \\
     \mathrm{exp} \left\{ -[(a - a_t)/a_c]^3 \right\} & a > a_t
     \end{array}
     \right.

For silicate grains :math:`D(a) = 0`.
For carbonaceous grains,

.. math::
  D(a) = \sum_{i=1}^2 \frac{B_i}{a} \mathrm{exp} \left\{ -\frac{1}{2}
     \left[ \frac{ln(a/a_{0,i})}{\sigma} \right]^2 \right\}

where

.. math::
  \begin{eqnarray}
  B_i & = & \frac{3}{(2\pi)^{3/2}} \frac{\mathrm{exp} (-4.5 \sigma^2)}{\rho a_{0,i}^3\sigma}
     \frac{b_{C,i} m_C}{1 + \mathrm{erf}(z)} \\
  z & = & \frac{3 \sigma}{\sqrt{2}} + \frac{\mathrm{ln} (a_{0,i}/a_{min})}{\sigma \sqrt{2}} \\
  \end{eqnarray}

and for carbonaceous material :math:`\sigma = 0.4`,
:math:`\rho = 2.24~\mathrm{cm}^{-3}`,
:math:`a_{0,1} = 0.35~\mathrm{nm}`, :math:`a_{0,2} = 3~\mathrm{nm}`,
:math:`a_{min} = 0.35~\mathrm{nm}`,
:math:`b_{C,1} = 0.75 b_C`, :math:`b_{C,2} = 0.25 b_C`, and
:math:`m_C` is the mass of a carbon atom.  Finally, :math:`b_C` is the
total C abundance in the two log-normal functions.

Example silicate and carbonaceous WD size distributions are shown below.

For the silicate grains, :math:`C = 1.33\times 10^{-11}`,
:math:`a_t = 171~\mathrm{nm}`, :math:`\alpha = -1.41`,
:math:`\beta = -11.5`, and :math:`a_c = 100~\mathrm{nm}`.

For the carbonaceous grains :math:`C = 4.15\times 10^{-11}`,
:math:`a_t = 8.37~\mathrm{nm}`, :math:`\alpha = -1.91`,
:math:`\beta = -0.125`, :math:`a_c = 499~\mathrm{nm}`,
:math:`b_C = 3\times 10^{-5}`.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from dgfit.dustmodel import WDDustModel

    fig, ax = plt.subplots()

    # simple intialization to get the WD size distribution function
    dmodel = WDDustModel()

    # size distributions in cm
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)*1e-8

    # silicate size distribution parameters
    cparams = [1.33e-12, 0.171e4, -1.41, -11.5]
    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, sizedist[indxs], label='silicate grains')

    # carbonaceous size distribution parameters
    cparams = [4.15e-11, 0.00837e4, -1.91, -0.125, 0.499e4, 3.0e-5]
    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, sizedist[indxs], label='carbonaceous grains')

    ax.set_ylim(1e-12, 100.)

    ax.set_xlabel('$a$ [$\mu m$]')
    ax.set_ylabel(r'$n_H^{-1} \; dn/da$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')

    plt.tight_layout()
    plt.show()

The same example size distributions are shown in a 2nd plot below
where the detailed wiggles in the overall powerlaw are emphasized
by multiplying by :math:`a^4`.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from dgfit.dustmodel import WDDustModel

    fig, ax = plt.subplots()

    # simple intialization to get the WD size distribution function
    dmodel = WDDustModel()

    # size distributions in cm
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)*1e-8

    # silicate size distribution parameters
    cparams = [1.33e-12, 0.171e4, -1.41, -11.5]
    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, 1e29*(a[indxs]**4)*sizedist[indxs], label='silicate grains')

    # carbonaceous size distribution parameters
    cparams = [4.15e-11, 0.00837e4, -1.91, -0.125, 0.499e4, 3.0e-5]
    sizedist = dmodel.compute_size_dist(a, cparams)

    # plot the nonzero sizedistribution values
    indxs, = np.where(sizedist > 0)
    ax.plot(a[indxs]*1e4, 1e29*(a[indxs]**4)*sizedist[indxs], label='carbonaceous grains')

    ax.set_ylim(0.01, 100.)

    ax.set_xlabel('$a$ [$\mu m$]')
    ax.set_ylabel(r'$10^{29} \; n_H^{-1} \; a^4 \; dn/da$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')

    plt.tight_layout()
    plt.show()
