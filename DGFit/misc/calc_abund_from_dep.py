#!/usr/bin/env python
# Started: Apr 2015 (KDG)
"""
Helper function to calculate abundances from depletion measurements
for individual lines-of-sight
"""
from __future__ import print_function, division

import numpy as np

from astropy.table import Table


def abundances_from_depletions(depletions, depletions_unc,
                               ref_abundances='SMC'):
    """
    abundances_from_depletions function

    Parameters
    ---------
    depletions: 'numpy.ndarray'
        depletions for C, O, Mg, Si, Fe
        units are log10[ N(X)/N(H)_gas / N(X)/N(H)_ref ]

    depletions_unc: 'numpy.ndarray'
        depletion uncertainties for C, O, Mg, Si, Fe
        units are log10[ N(X)/N(H)_gas / N(X)/N(H)_ref ]

    ref_abundances: 'string'
        string giving the source of the reference abundances
        default = 'SMC': allowed = ['LMC','SMC']

    Returns
    -------
    abundances: 'numpy.ndarray'
        abundances in dust grains for same elements
        units are N(X)/[10^6 N(H)]

    abundances_uncs: 'numpy.ndarray'
        abundance uncertainties

    ref_abund: 'numpy.ndarray'
        same but for the reference (gas+dust) abundances

    ref_abund_uncs: 'numpy.ndarray'
        reference (gas+dust) abundance uncertainties
    """
    if ref_abundances not in ['LMC', 'SMC']:
        print('ref_abundances of ' + ref_abundances + ' is not allowed')
        exit()

    # reference abundances
    ref_abund = {}
    ref_abund_unc = {}
    # need to add MW from Jenkins (2009)
    # LMC/SMC from Tchernyshyov et al. (2015?)
    ref_abund['LMC'] = np.array([7.94, 8.50, 7.26, 7.35, 7.32])
    ref_abund['SMC'] = np.array([7.52, 8.14, 6.88, 6.96, 6.89])
    ref_abund_unc['SMC'] = np.array([np.sqrt(0.10**2 + 0.04**2),
                                     np.sqrt(0.08**2 + 0.04**2),
                                     np.sqrt(0.06**2 + 0.03**2),
                                     np.sqrt(0.07**2 + 0.09**2),
                                     np.sqrt(0.08**2 + 0.03**2)])

    # compute the abundances
    # abund = (np.power(10.0,ref_abund[ref_abundances])*
    #         (1.0 - np.power(10.0,depletions)))

    # compute the uncertainties
    # C, O, Mg, Si, Fe
    tot_abund_log = ref_abund[ref_abundances]
    tot_abund_unc_log = ref_abund_unc[ref_abundances]

    tot_abund = np.power(10., tot_abund_log - 6.0)
    tot_abund_up = np.power(10., tot_abund_log + tot_abund_unc_log - 6.0)
    tot_abund_down = np.power(10., tot_abund_log - tot_abund_unc_log - 6.0)
    tot_abund_unc = 0.5*(tot_abund_up - tot_abund_down)

    gas_abund_log = depletions
    gas_abund_unc_log = depletions_unc

    gas_abund = tot_abund*np.power(10.0, gas_abund_log)
    gas_abund_up = tot_abund*np.power(10.0, gas_abund_log + gas_abund_unc_log)
    gas_abund_down = tot_abund*np.power(10.0,
                                        gas_abund_log - gas_abund_unc_log)
    gas_abund_unc = 0.5*(gas_abund_up - gas_abund_down)

    dust_abund = tot_abund - gas_abund
    dust_abund_unc = (np.sqrt(np.square(tot_abund_unc)
                              + np.square(gas_abund_unc)))

    return (dust_abund, dust_abund_unc, tot_abund, tot_abund_unc)


if __name__ == "__main__":

    atomnames = ['C', 'O', 'Mg', 'Si', 'Fe']
    dep_smc_azv215 = np.array([-0.23, -0.95, -0.89, -0.52, -1.87])
    dep_unc_smc_azv215 = np.array([0.15, 0.60, 0.15, 0.76, 0.11])

    results = abundances_from_depletions(dep_smc_azv215, dep_unc_smc_azv215,
                                         ref_abundances='SMC')

    # number of atoms per 10^6 H atoms
    print(atomnames)
    print(results[0])
    print(results[1])
    print(results[2])
    print(results[3])

    # write a table in the format desired by DGFit

    a = Table([atomnames, results[0], results[1], results[2], results[3]],
              names=('atom', 'abund', 'abund_unc', 'total_abund',
                     'total_abund_unc'))

    a.write('SMC_AzV215.dat', format='ascii.commented_header')
