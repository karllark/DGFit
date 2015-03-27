#!/usr/bin/env python
#
# dg_props objects
#  dust grain properties stored by dust size/composition
#
# Started: Jan 2015 (KDG)
# 
from __future__ import print_function
import glob
import sys
import re
import string
import math

#import numexpr

import numpy as np
from astropy.table import Table

import matplotlib.pyplot as pyplot

# Object for the proprerties of dust grain with a specific composition
class DustGrains():
    def __init__(self, componentname, path='./',
                 min_wave=0., max_wave=1e6,
                 min_wave_emission=0., max_wave_emission=1e6):
        # check that the component name is allowed
        #_allowed_components = ['astro-silicates','astro-graphite','astro-carbonaceous','astro-PAH']
        _allowed_components = ['astro-silicates','astro-carbonaceous']
        if not componentname in _allowed_components:
            print(componentname + ' not one of the allowed grain components')
            print(allowed_components)
            exit()

        # set useful quantities for each composition
        if componentname == 'astro-silicates':   # from WD01
            self.density = 3.5  # g/cm^3
            self.atomic_composition = 'MgFeSiO4'
            self.atomic_comp_names = ['Mg','Fe','Si','O']
            self.atomic_comp_number = np.array([1, 1, 1, 4])
            self.atomic_comp_masses = np.array([24.305, 55.845, 28.0855, 15.994])*1.660e-24  # in grams
        elif componentname == 'astro-carbonaceous':  # from WD01
            self.density = 2.24  # g/cm^3  
            self.atomic_composition = 'C'
            self.atomic_comp_names = ['C']
            self.atomic_comp_number = np.array([1])
            self.atomic_comp_masses = np.array([12.0107])*1.660e-24  # in grams

        #useful quantities
        self.mass_per_mol_comp = np.sum(self.atomic_comp_masses*self.atomic_comp_number)
        self.col_den_constant = (4./3.)*math.pi*self.density*self.atomic_comp_number/self.mass_per_mol_comp

        # get the filenames of this component for all sizes
        filelist = []
        sizenum = -1
        for file in glob.glob(path+"INDIV-GRAINS-FAKE-FIT_c_*"+componentname+"*.dat"):
            m = re.search('_s_(.+?).dat', file)
            if m:
                found = m.group(1)
                sizenum = found

                # get the grain size
                f = open(file, 'r')
                firstline = f.readline()
                space_pos = string.find(firstline,' ',5)
                secondline = f.readline()
                colon_pos = string.find(secondline, ':')
                f.close()

                if secondline[colon_pos+2:colon_pos+5] == 'Yes':
                    stochastic_heating = True
                else:
                    stochastic_heating = False
                
                filelist.append((file,int(sizenum),float(firstline[1:space_pos]),stochastic_heating))

        # temp code to just pick every 10th size
        tindxs = np.arange(0,len(filelist),10)
        sfilelist = sorted(filelist, key=lambda file: file[1])
        filelist = []
        for k in tindxs:
            filelist.append(sfilelist[k])

        # setup the variables to store the grain information
        self.name = componentname
        self.n_sizes = len(filelist)
        self.sizes = np.empty(self.n_sizes)
        self.size_dist = np.empty(self.n_sizes)
        self.stochastic_heating = np.empty(self.n_sizes)

        # loop over the files from the smallest to the largest sizes
        for k, file in enumerate(sorted(filelist, key=lambda file: file[1])):
            # read in the table of grain properties for this size
            t = Table.read(file[0],format='ascii.commented_header',header_start=9)

            # setup more variables now that we know the number of wavelengths
            if k == 0:
                # generate the indices to crop the wavelength to the desired range
                gindxs, = np.where((t['Wavelength'] >= min_wave) & (t['Wavelength'] <= max_wave))
                egindxs, = np.where((t['Wavelength'] >= min_wave_emission) & (t['Wavelength'] <= max_wave_emission))
                self.wavelengths = t['Wavelength'][gindxs]
                self.wavelengths_emission = t['Wavelength'][egindxs]
                self.n_wavelengths = len(self.wavelengths)
                self.n_wavelengths_emission = len(self.wavelengths_emission)
                self.cext = np.empty((self.n_sizes, self.n_wavelengths))
                self.cabs = np.empty((self.n_sizes, self.n_wavelengths))
                self.csca = np.empty((self.n_sizes, self.n_wavelengths))
                self.emission = np.empty((self.n_sizes, self.n_wavelengths_emission))

            # store the info
            self.sizes[k] = file[2]
            self.stochastic_heating[k] = file[3]
            self.cext[k,:] = t['CExt'][gindxs]
            self.csca[k,:] = t['CSca'][gindxs]
            self.cabs[k,:] = t['CAbs'][gindxs]
            if file[3]:
                self.emission[k,:] = t['StEmission'][egindxs]
            else:
                self.emission[k,:] = t['EqEmission'][egindxs]

            # convert emission from ergs/(s cm sr) to Jy/sr
            #   wavelengths in microns
            self.emission[k,:] *= (self.wavelengths_emission)**2/2.998e10  # convert from cm^-1 to Hz^-1
            self.emission[k,:] /= 1e-19  # convert from ergs/(s Hz) to Jy
            self.emission[k,:] *= 1e-6 # convert from Jy/sr to MJy/sr
            self.emission[k,:] *= 1e-4 # arbitrary value needed to match observations - need to find the origin
            
            # default size distributions - MRN distribution
            self.size_dist[k] = self.sizes[k]**(-3.5)

    # function to integrate this component
    # returns the effective/total cabs, csca, etc.
    # these are normalized to NHI (assumption)
    def eff_grain_props(self):
        # initialize the results
        _effcabs = np.empty(self.n_wavelengths)
        _effcsca = np.empty(self.n_wavelengths)
        _natoms = np.empty(len(self.atomic_comp_names))
        _emission = np.empty(self.n_wavelengths_emission)
        # do a very simple integration (later this could be made more complex)
        deltas = 0.5*(self.sizes[1:self.n_sizes] - self.sizes[0:self.n_sizes-1])
        sizedist1 = self.size_dist[0:self.n_sizes-1]
        sizedist2 = self.size_dist[1:self.n_sizes]
        for i in range(self.n_wavelengths):
            _effcabs[i] = np.sum( deltas*( (self.cabs[0:self.n_sizes-1,i]*sizedist1) + (self.cabs[1:self.n_sizes,i]*sizedist2) ) )
            _effcsca[i] = np.sum( deltas*( (self.csca[0:self.n_sizes-1,i]*sizedist1) + (self.csca[1:self.n_sizes,i]*sizedist2) ) )

            # *not* faster
            #_effcabs[i] = np.sum( numexpr.evaluate('deltas*( (cabs1*sizedist1) + (cabs2*sizedist2) )',
            #                                       local_dict={'deltas':deltas,
            #                                                   'cabs1':self.cabs[0:self.n_sizes-1,i],
            #                                                   'sizedist1':self.size_dist[0:self.n_sizes-1],
            #                                                   'cabs2':self.cabs[1:self.n_sizes,i],
            #                                                   'sizedist2':self.size_dist[1:self.n_sizes]} ) )

            #for k in range(self.n_sizes-1):
            #    delta = 0.5*(self.sizes[k+1] - self.sizes[k])  # factor of 0.5 is here for speed

        # compute the integrated emission spectrum
        for i in range(self.n_wavelengths_emission):
            _emission[i] = np.sum( deltas*( (self.emission[0:self.n_sizes-1,i]*sizedist1) + (self.emission[1:self.n_sizes,i]*sizedist2) ) )

        # compute the number of atoms/NHI 
        for i in range(len(self.atomic_comp_names)):
            _natoms[i] = np.sum(deltas*( ((self.sizes[0:self.n_sizes-1]**3)*self.size_dist[0:self.n_sizes-1]*self.col_den_constant[i]) +
                                         ((self.sizes[1:self.n_sizes]**3)*self.size_dist[1:self.n_sizes]*self.col_den_constant[i]) ) )

        #for k in range(self.n_sizes-1):
        #    delta = 0.5*(self.sizes[k+1] - self.sizes[k])  # factor of 0.5 is here for speed
        #    for i in range(len(self.atomic_comp_names)):
        #        _natoms[i] += delta*( ((self.sizes[k+1]**3)*self.size_dist[k+1]*self.col_den_constant[i]) +
        #                              ((self.sizes[k]**3)*self.size_dist[k]*self.col_den_constant[i]) )

        # convert to N(N) per 1e6 N(HI)
        _natoms *= 1e6
            
        # return the results as a tuple of arrays
        return (_effcabs, _effcsca, dict(zip(self.atomic_comp_names, _natoms)), _emission)

if __name__ == "__main__":
    
    dustgraincomp = DustGrains('astro-silicates',path='/home/kgordon/Dirty_v2/write_grain/indiv_grain/',
                               min_wave=0.09,max_wave=3.0)

    cabs, csca, natoms = dustgraincomp.eff_grain_props()

    print(natoms)

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(dustgraincomp.wavelengths, cabs+csca, 'g-')
    ax.plot(dustgraincomp.wavelengths, cabs, 'r-')
    ax.plot(dustgraincomp.wavelengths, csca, 'b-')
    ax.plot(dustgraincomp.wavelengths, csca/(cabs+csca), 'c-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    pyplot.show()
    exit()
