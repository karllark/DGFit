#!/usr/bin/env python
#
# Observed Data Object
#
# Started: Jan 2015 (KDG)
# 
from __future__ import print_function

import string

import numpy as np
import matplotlib.pyplot as pyplot
from astropy.table import Table

from Ext import f99
from Ext import extdata

# Object for the observed dust data
class ObsData():
    def __init__(self, wavelengths, wavelengths_emission):
        # read in the observed extinction curve
        extinfo = extdata.ExtData()
        extinfo.read_ext_data('Ext/azv18_ext.fits')

        # compute the dust extintion curve one the supplied wavelength grid
        self.n_wavelengths = len(wavelengths)
        self.wavelengths = wavelengths
        self.Rv = extinfo.rv[0]
        self.avnhi = 1./extinfo.nhiav[0]
        self.alnhi = self.avnhi*f99.f99(self.Rv, 1.0/wavelengths, x0=extinfo.fm90['x0'][0], gamma=extinfo.fm90['gamma'][0],
                                        c2=extinfo.fm90['C2'][0],c3=extinfo.fm90['C3'][0],c4=extinfo.fm90['C4'][0])

        # Kirill sent those column densities toward AzV 18:
        # Fe II : 15.58 +/- 0.02
        # Mg II : 16.50 +/- 0.03
        # O I : 18.27 +/- 0.09
        # Si II : 16.63 +/- 0.03
        #
        # Kirill sent at depletion for carbon of -0.15  (F_* = 0.5)

        # C, O, Mg, Si, Fe
        tot_abund_log = np.array([7.52, 8.14, 6.88, 6.96, 6.89])
        tot_abund_unc_log = np.array([0.14, 0.12, 0.09, 0.16, 0.11])

        tot_abund = np.power(10.,tot_abund_log - 6.0)
        tot_abund_up = np.power(10.,tot_abund_log + tot_abund_unc_log - 6.0)
        tot_abund_down = np.power(10.,tot_abund_log - tot_abund_unc_log - 6.0)
        tot_abund_unc = 0.5*(tot_abund_up - tot_abund_down)
        #print(tot_abund, tot_abund_unc)

        gas_abund_log = np.array([18.27, 16.50, 16.63, 15.58])
        gas_abund_unc_log = np.array([0.09, 0.03, 0.03, 0.02])
        gas_abund = 1e6*np.power(10.,gas_abund_log)/np.power(10.,22.04)
        nhi = np.power(10.,22.04)
        gas_abund_up = 1e6*np.power(10.,gas_abund_log+gas_abund_unc_log)/nhi
        gas_abund_down = 1e6*np.power(10.,gas_abund_log-gas_abund_unc_log)/nhi
        gas_abund_unc = 0.5*(gas_abund_up - gas_abund_down)
        #print(gas_abund_unc)

        nhi_unc = 0.5*(np.power(10.,22.04+0.12) - np.power(10.,22.04-0.18))

        gas_abund = np.concatenate(([1e-6*np.power(10.,tot_abund_log[0]-0.15)],gas_abund))
        gas_abund_unc = np.concatenate(([0.2*gas_abund[0]],gas_abund_unc))

        gas_abund_unc = gas_abund*np.sqrt((np.power((gas_abund_unc/gas_abund),2) + np.power((nhi_unc/nhi),2)))

        #print(gas_abund, gas_abund_unc)

        dust_abund = tot_abund - gas_abund
        dust_abund_unc = np.sqrt(np.square(tot_abund_unc) + np.square(gas_abund_unc))

        #print(dust_abund, dust_abund_unc)

        self.total_abundance = {'C': (tot_abund[0], tot_abund_unc[0]),
                                'O': (tot_abund[1], tot_abund_unc[1]),
                                'Mg': (tot_abund[2], tot_abund_unc[2]),
                                'Si': (tot_abund[3], tot_abund_unc[3]),
                                'Fe': (tot_abund[4], tot_abund_unc[4])}

        self.depletions = {'C': (dust_abund[0], dust_abund_unc[0]),
                           'O': (dust_abund[1], dust_abund_unc[1]),
                           'Mg': (dust_abund[2], dust_abund_unc[2]),
                           'Si': (dust_abund[3], dust_abund_unc[3]),
                           'Fe': (dust_abund[4], dust_abund_unc[4])}

        self.fit_depletions = True
        self.fit_ir_emission = False
        self.fit_scat_param = False

if __name__ == "__main__":
    
    print('No test code')
