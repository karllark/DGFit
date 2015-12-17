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

# Object for the observed dust data
class ObsData():
    def __init__(self, wavelengths, wavelengths_emission):
        # compute the dust extintion curve one the supplied wavelength grid
        self.n_wavelengths = len(wavelengths)
        self.wavelengths = wavelengths
        self.Rv = 3.1
        self.avnhi = self.Rv/5.8e21   # number is N(HI)/E(B-V) from Bohlin, Savage, & Drake (1978)
        self.alnhi = self.avnhi*f99.f99(self.Rv, 1.0/wavelengths)

        # depletions from Jenkins (2009) for F* = 0.36 in N(X) per 1e6 N(HI)
        self.depletions = {'C': (83.41, 0.05*83.41), 'O': (109.26, 0.05*109.26), 'Mg': (31.90, 0.05*31.90),
                           'Si': (31.24, 0.05*31.24), 'Fe': (33.34, 0.05*33.34)}

        self.total_abundance = {'C': (2*83.41, 0.05*2*83.41), 'O': (2*109.26, 0.05*2*109.26), 'Mg': (2*31.90, 0.05*2*31.90),
                                'Si': (2*31.24, 0.05*2*31.24), 'Fe': (2*33.34, 0.05*2*33.34)}

        # diffuse IR emission spectrum (Gordon et al. 2014)
        self.ir_emission_waves = np.array([100., 160., 250., 350., 500.])
        self.ir_emission = np.array([0.71, 1.53, 1.08, 0.56, 0.25])/1e20
        self.ir_emission_uncs = np.array([0.71, 1.53, 1.08, 0.56, 0.25])*0.05/1e20
        # get the indxs to the nearest wavelenths in the model emission spectra
        self.n_ir_emission_waves = len(self.ir_emission_waves)
        self.ir_emission_indxs = np.empty(self.n_ir_emission_waves, np.int)
        for i in range(self.n_ir_emission_waves):
            sindxs = np.argsort(np.abs(self.ir_emission_waves[i] - wavelengths_emission))
            self.ir_emission_indxs[i] = sindxs[0]

        # dust albedo (Gordon et al. AoD proceedings)
        files_dgl = ["mathis73","morgan76","lillie76","toller81", 
                     "murthy93","murthy95","petersohn97","witt97", 
                     "schiminovich01","shalima04","sujatha05","sujatha07","sujatha10"]
        scat_path = "/home/kgordon/Pro/Dust/Scat_Data/"

        scat_waves = []
        scat_albedo = []
        scat_albedo_unc = []
        scat_g = []
        scat_g_unc = []
        scat_ref = []
        scat_indxs = []
        for sfile in files_dgl:
            f = open(scat_path + sfile + '.dat', 'r')
            ref = f.readline().rstrip()
            f.close()

            t = Table.read(scat_path+sfile+'.dat',format='ascii',header_start=1)
            for k in range(len(t)):
                scat_waves.append(t['wave,'][k])
                scat_albedo.append(t['albedo,'][k])
                scat_albedo_unc.append(t['delta,'][k])
                scat_g.append(t['g,'][k])
                scat_g_unc.append(t['delta'][k])
                scat_ref.append(ref)
                # get the index of the nearest wavelength in the model grid
                sindxs = np.argsort(np.abs(scat_waves[-1] - wavelengths))
                scat_indxs.append(sindxs[0])

        self.scat_waves = np.array(scat_waves)*1e-4
        self.scat_albedo = np.array(scat_albedo)
        self.scat_albedo_unc = np.array(scat_albedo_unc)
        self.scat_g = np.array(scat_g)
        self.scat_g_unc = np.array(scat_g_unc)
        self.scat_ref = scat_ref
        self.scat_indxs = np.array(scat_indxs, np.int)

        # get the indxs to the nearest wavelenths in the model emission spectra
        self.n_ir_emission_waves = len(self.ir_emission_waves)
        self.ir_emission_indxs = np.empty(self.n_ir_emission_waves, np.int)
        for i in range(self.n_ir_emission_waves):
            sindxs = np.argsort(np.abs(self.ir_emission_waves[i] - wavelengths_emission))
            self.ir_emission_indxs[i] = sindxs[0]

        self.fit_depletions = True
        self.fit_ir_emission = True
        self.fit_scat_param = True
        
if __name__ == "__main__":
    
    print('No test code')
