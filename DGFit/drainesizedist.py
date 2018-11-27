#!/usr/bin/env python
#
# generate a Draine dust grain size distribution
#   from Weingartner & Draine (2001)
#
# Started: Jul 2016 (KDG)

from scipy.special import erf
import numpy as np

from DGFit.DustModel import WDDustModel


def WGsizedist_graphite(a, bC_input, Cg, a_tg, a_cg, alpha_g, beta_g):
    # a in Angstrom and assumed to always be > 3.5 A

    # very small gain size distribution
    a0 = np.array([3.5, 30.])   # in A
    bC = np.array([0.75, 0.25])*bC_input
    sigma = 0.4
    rho = 2.24  # in g/cm^3 for graphite
    mC = 12.0107*1.660e-24

    Da = 0.0
    for i in range(2):
        Bi = ((3.0/(np.power(2.0*np.pi, 1.5)))
              * (np.exp(-4.5*np.power(sigma, 2.0))
              / (rho*np.power(1e-8*a0[i], 3.0)*sigma))
              * (bC[i]*mC/(1.0 + erf((3.0*sigma/np.sqrt(2.0))
                                     + np.log(a0[i]/3.5)/(sigma*np.sqrt(2.0)))
                           )))

        Da += (Bi/(1e-8*a))*np.exp(-0.5*np.power(np.log(a/a0[i])/sigma, 2.0))

    # larger grain size distribution
    if beta_g >= 0.0:
        Fa = 1.0 + beta_g*a/a_tg
    else:
        Fa = 1.0/(1.0 - beta_g*a/a_tg)

    Ga = np.full((len(a)), 1.0)
    indxs, = np.where(a > a_tg)
    Ga[indxs] = np.exp(-1.0*np.power((a[indxs] - a_tg)/a_cg, 3.0))

    graphite_sd = Da + (Cg/(1e-8*a))*np.power(a/a_tg, alpha_g)*Fa*Ga

    return graphite_sd


if __name__ == "__main__":

    import argparse

    import matplotlib.pyplot as plt
    import matplotlib

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    # grain sizes
    a = np.logspace(np.log10(3.5), np.log10(2e4), num=200)

    # setup the plots
    fontsize = 18
    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = plt.subplots(ncols=2, figsize=(14, 8))

    # coefficients from Weingartner & Draine for R(V) = 3.1, Case A
    bC = np.array([0.0, 1.0, 2.0, 3.0,
                   4.0, 5.0, 6.0])*1e-5
    Cg = np.array([9.94e-11, 3.79e-10, 5.57e-11, 4.15e-11,
                   2.90e-11, 3.20e-12, 9.99e-12])
    atg = np.array([0.00745, 0.00373, 0.00828, 0.00837,
                    0.00898, 0.0254, 0.0107])*1e4
    acg = np.array([0.606, 0.586, 0.543, 0.499,
                    0.489, 0.438, 0.428])*1e4
    alpha_g = np.array([-2.25, -2.17, -2.04, -1.91,
                        -1.84, -1.72, -1.54])
    beta_g = np.array([-0.0648, -0.0382, -0.111, -0.125,
                       -0.132, -0.322, -0.165])

    dmodel = WDDustModel()

    # get the size distributions and plot
    for i in range(len(bC)):
        params = [Cg[i], atg[i], alpha_g[i], beta_g[i], acg[i], bC[i]]
        print(params)
        dm_sizedist_graphite = dmodel.compute_size_dist(a*1e-8, params)
        sizedist_graphite = WGsizedist_graphite(a,
                                                bC[i],
                                                Cg[i],
                                                atg[i],
                                                acg[i],
                                                alpha_g[i],
                                                beta_g[i])

        ax[0].plot(a*1e-4, sizedist_graphite*np.power(1e-8*a, 4.0)*1e29)
        ax[1].plot(a*1e-4, sizedist_graphite)
        ax[1].plot(a*1e-4, dm_sizedist_graphite, 'k--')

    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlim(3.5e-4, 2.2)
    ax[0].set_ylim(0.1, 25)
    ax[0].set_xlabel(r'$a$ [$\mu m$]')
    ax[0].set_ylabel(r'$dn/da$ [$10^{-29}$ $cm^3$ $n_H^{-1}$]')

    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlim(3.5e-4, 2.2)
    ax[1].set_ylim(1e-15, 1e4)
    ax[1].set_xlabel(r'$a$ [$\mu m$]')
    ax[1].set_ylabel(r'$dn/da$ [$n_H^{-1}$]')
    plt.tight_layout()

    # show or save
    basename = 'drainesizedist'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        plt.show()
