#!/usr/bin/env python
#
# dg_props objects
#  dust grain properties stored by dust size/composition
#
# Started: Jan 2015 (KDG)
# Updated to include better diagnoistic plots when run (Mar 2016 KDG)
# 
from __future__ import print_function
import glob
import sys
import re
import string
import math
import colorsys
import argparse

import numpy as np
from astropy.table import Table

import matplotlib.pyplot as pyplot
import matplotlib

from scipy.interpolate import interp1d

from ObsData import ObsData

# Object for the proprerties of dust grain with a specific composition
class DustGrains():
    def __init__(self):
        self.origin = None

    def from_files(self, componentname, path='./',
                   min_wave=0., max_wave=1e6,
                   min_wave_emission=0., max_wave_emission=1e6):
        self.origin = 'files'
        
        # check that the component name is allowed
        #_allowed_components = ['astro-silicates','astro-graphite',
        #                       'astro-carbonaceous','astro-PAH']
        _allowed_components = ['astro-silicates','astro-carbonaceous',
                               'astro-graphite']
        if not componentname in _allowed_components:
            print(componentname + ' not one of the allowed grain components')
            print(_allowed_components)
            exit()

        # set useful quantities for each composition
        if componentname == 'astro-silicates':   # from WD01
            self.density = 3.5  # g/cm^3
            self.atomic_composition = 'MgFeSiO4'
            self.atomic_comp_names = ['Mg','Fe','Si','O']
            self.atomic_comp_number = np.array([1, 1, 1, 4])
            self.atomic_comp_masses = np.array([24.305, 55.845, 28.0855,
                                                15.994])*1.660e-24  # in grams
        elif componentname == 'astro-carbonaceous':  # from WD01
            self.density = 2.24  # g/cm^3  
            self.atomic_composition = 'C'
            self.atomic_comp_names = ['C']
            self.atomic_comp_number = np.array([1])
            self.atomic_comp_masses = np.array([12.0107])*1.660e-24  # in grams

        elif componentname == 'astro-graphite':  # need origin (copy)
            self.density = 2.24  # g/cm^3  
            self.atomic_composition = 'C'
            self.atomic_comp_names = ['C']
            self.atomic_comp_number = np.array([1])
            self.atomic_comp_masses = np.array([12.0107])*1.660e-24  # in grams
            
        #useful quantities
        self.mass_per_mol_comp = np.sum(self.atomic_comp_masses*
                                        self.atomic_comp_number)
        self.col_den_constant = (4./3.)*math.pi*self.density* \
                                self.atomic_comp_number/self.mass_per_mol_comp

        # get the filenames of this component for all sizes
        filelist = []
        sizenum = -1
        for file in glob.glob(path+"INDIV-GRAINS-FAKE-FIT_c_*"+
                              componentname+"*.dat"):
            m = re.search('_s_(.+?).dat', file)
            if m:
                found = m.group(1)
                sizenum = found

                # get the grain size
                f = open(file, 'r')
                firstline = f.readline()
                space_pos = firstline.find(' ',5)
                secondline = f.readline()
                colon_pos = secondline.find( ':')
                f.close()

                if secondline[colon_pos+2:colon_pos+5] == 'Yes':
                    stochastic_heating = True
                else:
                    stochastic_heating = False
                
                filelist.append((file,int(sizenum),
                                 float(firstline[1:space_pos]),
                                 stochastic_heating))

        # temp code to just pick every 5th size
        tindxs = np.arange(0,len(filelist),5)
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
            t = Table.read(file[0],format='ascii.commented_header',
                           header_start=9)
            
            # setup more variables now that we know the number of wavelengths
            if k == 0:
                # generate the indices to crop the wavelength to the
                #      desired range
                gindxs, = np.where((t['Wavelength'] >= min_wave) &
                                   (t['Wavelength'] <= max_wave))
                egindxs, = np.where((t['Wavelength'] >= min_wave_emission) &
                                    (t['Wavelength'] <= max_wave_emission))
                self.wavelengths = np.array(t['Wavelength'][gindxs])
                self.wavelengths_emission = np.array(t['Wavelength'][egindxs])
                self.n_wavelengths = len(self.wavelengths)
                self.n_wavelengths_emission = len(self.wavelengths_emission)
                self.cext = np.empty((self.n_sizes, self.n_wavelengths))
                self.cabs = np.empty((self.n_sizes, self.n_wavelengths))
                self.csca = np.empty((self.n_sizes, self.n_wavelengths))
                self.scat_g = np.empty((self.n_sizes, self.n_wavelengths))
                self.emission = np.empty((self.n_sizes,
                                          self.n_wavelengths_emission))

            # store the info
            self.sizes[k] = file[2]
            self.stochastic_heating[k] = file[3]
            self.cext[k,:] = t['CExt'][gindxs]
            self.csca[k,:] = t['CSca'][gindxs]
            self.cabs[k,:] = t['CAbs'][gindxs]
            self.scat_g[k,:] = t['G'][gindxs]
            if file[3]:
                self.emission[k,:] = t['StEmission'][egindxs]
            else:
                self.emission[k,:] = t['EqEmission'][egindxs]

            # convert emission from ergs/(s cm sr) to Jy/sr
            #   wavelengths in microns
            #      convert from cm^-1 to Hz^-1
            self.emission[k,:] *= (self.wavelengths_emission)**2/2.998e10  
            self.emission[k,:] /= 1e-19  # convert from ergs/(s Hz) to Jy
            self.emission[k,:] *= 1e-6 # convert from Jy/sr to MJy/sr
            # convert from m^-2 to cm^-2
            self.emission[k,:] *= 1e-4 
            
            # default size distributions - MRN distribution
            self.size_dist[k] = self.sizes[k]**(-3.5)

        # aliases for albedo and g calculations
        #    here they are on the same wavelength grid
        #    when calculated from an ObsData object they are not
        #    see the next function
        #    (done to allow rest of DustGrain code to be generic)
        self.n_wavelengths_scat_a = self.n_wavelengths
        self.wavelengths_scat_a = self.wavelengths
        self.scat_a_cext = self.cext
        self.scat_a_csca = self.csca

        self.n_wavelengths_scat_g = self.n_wavelengths
        self.wavelengths_scat_g = self.wavelengths
        self.scat_g_csca = self.csca

    # generate the dust grain info on the observed data grids
    #   based on an existing DustGrain object
    def from_object(self, DustGrain, ObsData):
        self.origin = 'object'

        # copy the basic information on the grain
        self.density = DustGrain.density
        self.atomic_composition = DustGrain.atomic_composition
        self.atomic_comp_names = DustGrain.atomic_comp_names
        self.atomic_comp_number = DustGrain.atomic_comp_number
        self.atomic_comp_masses = DustGrain.atomic_comp_masses
        self.mass_per_mol_comp = DustGrain.mass_per_mol_comp
        self.col_den_constant = DustGrain.col_den_constant

        self.name = DustGrain.name
        self.n_sizes = DustGrain.n_sizes
        self.sizes = DustGrain.sizes
        self.size_dist = DustGrain.size_dist
        self.stochastic_heating = DustGrain.stochastic_heating

        # new values on the observed wavelength grids
        self.wavelengths = ObsData.ext_waves
        self.wavelengths_emission = ObsData.ir_emission_waves
        self.wavelengths_scat_a = ObsData.scat_a_waves
        self.wavelengths_scat_g = ObsData.scat_g_waves
        self.n_wavelengths = len(self.wavelengths)
        self.n_wavelengths_emission = len(self.wavelengths_emission)
        self.n_wavelengths_scat_a = len(self.wavelengths_scat_a)
        self.n_wavelengths_scat_g = len(self.wavelengths_scat_g)

        # variables to store the dust grain properties
        self.cext = np.empty((self.n_sizes, self.n_wavelengths))
        self.cabs = np.empty((self.n_sizes, self.n_wavelengths))
        self.csca = np.empty((self.n_sizes, self.n_wavelengths))
        self.scat_a_cext = np.empty((self.n_sizes, self.n_wavelengths_scat_a))
        self.scat_a_csca = np.empty((self.n_sizes, self.n_wavelengths_scat_a))
        self.scat_g = np.empty((self.n_sizes, self.n_wavelengths_scat_g))
        self.scat_g_csca = np.empty((self.n_sizes, self.n_wavelengths_scat_g))
        self.emission = np.empty((self.n_sizes,
                                  self.n_wavelengths_emission))

        # loop over the sizes and generate grain info on the observed data grid
        for i in range(self.n_sizes):
            cext_interp = interp1d(DustGrain.wavelengths, DustGrain.cext[i,:])
            cabs_interp = interp1d(DustGrain.wavelengths, DustGrain.cabs[i,:])
            csca_interp = interp1d(DustGrain.wavelengths, DustGrain.csca[i,:])
            self.cext[i,:] = cext_interp(self.wavelengths)
            self.cabs[i,:] = cabs_interp(self.wavelengths)
            self.csca[i,:] = csca_interp(self.wavelengths)

            self.scat_a_cext[i,:] = cext_interp(self.wavelengths_scat_a)
            self.scat_a_csca[i,:] = csca_interp(self.wavelengths_scat_a)

            g_interp = interp1d(DustGrain.wavelengths, DustGrain.scat_g[i,:])
            self.scat_g[i,:] = g_interp(self.wavelengths_scat_g)
            self.scat_g_csca[i,:] = csca_interp(self.wavelengths_scat_g)

            emission_interp = interp1d(DustGrain.wavelengths_emission, 
                                       DustGrain.emission[i,:])
            self.emission[i,:] = emission_interp(self.wavelengths_emission)


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
        deltas = 0.5*(self.sizes[1:self.n_sizes] -
                      self.sizes[0:self.n_sizes-1])
        sizedist1 = self.size_dist[0:self.n_sizes-1]
        sizedist2 = self.size_dist[1:self.n_sizes]
        for i in range(self.n_wavelengths):
            _effcabs[i] = np.sum( deltas*( (self.cabs[0:self.n_sizes-1,i]*
                                            sizedist1) +
                                           (self.cabs[1:self.n_sizes,i]*
                                            sizedist2) ) )
            _effcsca[i] = np.sum( deltas*( (self.csca[0:self.n_sizes-1,i]*
                                            sizedist1) +
                                           (self.csca[1:self.n_sizes,i]*
                                            sizedist2) ) )

            # *not* faster to use numexpr (tested in 2015)

        # scattering parameters a & g
        n_waves_scat_a = self.n_wavelengths_scat_a
        scat_a_cext = self.scat_a_cext
        scat_a_csca = self.scat_a_csca

        _effscat_a_cext = np.empty(n_waves_scat_a)
        _effscat_a_csca = np.empty(n_waves_scat_a)

        for i in range(n_waves_scat_a):
            _effscat_a_cext[i] = np.sum( deltas*( 
                    (scat_a_cext[0:self.n_sizes-1,i]*
                     sizedist1) +
                    (scat_a_cext[1:self.n_sizes,i]*
                     sizedist2) ) )
            _effscat_a_csca[i] = np.sum( deltas*( 
                    (scat_a_csca[0:self.n_sizes-1,i]*
                     sizedist1) +
                    (scat_a_csca[1:self.n_sizes,i]*
                     sizedist2) ) )

        n_waves_scat_g = self.n_wavelengths_scat_a
        scat_g_csca = self.scat_g_csca

        _effg = np.empty(n_waves_scat_g)
        _effscat_g_csca = np.empty(n_waves_scat_g)

        for i in range(n_waves_scat_g):
            _effg[i] = np.sum( deltas*( 
                    (self.scat_g[0:self.n_sizes-1,i]*
                     scat_g_csca[0:self.n_sizes-1,i]*
                     sizedist1) +
                    (self.scat_g[1:self.n_sizes,i]*
                     scat_g_csca[1:self.n_sizes,i]*
                     sizedist2) ) )
            _effscat_g_csca[i] = np.sum( deltas*( 
                    (scat_g_csca[0:self.n_sizes-1,i]*
                     sizedist1) +
                    (scat_g_csca[1:self.n_sizes,i]*
                     sizedist2) ) )

        # compute the integrated emission spectrum
        for i in range(self.n_wavelengths_emission):
            _emission[i] = np.sum( deltas*( (self.emission[0:self.n_sizes-1,i]*
                                             sizedist1) +
                                            (self.emission[1:self.n_sizes,i]*
                                             sizedist2) ) )

        # compute the number of atoms/NHI 
        for i in range(len(self.atomic_comp_names)):
            _natoms[i] = np.sum(deltas*( ((self.sizes[0:self.n_sizes-1]**3)*
                                          self.size_dist[0:self.n_sizes-1]*
                                          self.col_den_constant[i]) +
                                         ((self.sizes[1:self.n_sizes]**3)*
                                          self.size_dist[1:self.n_sizes]*
                                          self.col_den_constant[i]) ) )

        # convert to N(N) per 1e6 N(HI)
        _natoms *= 1e6
            
        # return the results as a tuple of arrays
        return (_effcabs,
                _effcsca,
                dict(zip(self.atomic_comp_names, _natoms)),
                _emission,
                _effscat_a_csca/_effscat_a_cext,
                _effg/_effscat_g_csca,
                _effscat_a_cext,
                _effscat_a_csca,
                _effscat_g_csca)

if __name__ == "__main__":
    
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition",
                        choices=['astro-silicates', 'astro-carbonaceous',
                                 'astro-graphite',
                                 'astro-PAH-ionized', 'astro-PAH-neutral'],
                        default='astro-silicates', 
                        help="Grain composition")
    parser.add_argument("--obsdata", help="transform to observed data grids",
                        action="store_true")
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    DG = DustGrains()
    DG.from_files(args.composition,
                  path='/home/kgordon/Dirty_v2/write_grain/indiv_grain2/')

    if args.obsdata:
        OD = ObsData(['data_mw_rv31/MW_diffuse_Gordon09_band_ext.dat',
                      'data_mw_rv31/MW_diffuse_Gordon09_iue_ext.dat',
                      'data_mw_rv31/MW_diffuse_Gordon09_fuse_ext.dat'],
                     'data_mw_rv31/MW_diffuse_Jenkins09_abundances.dat',
                     'data_mw_rv31/MW_diffuse_Compiegne11_ir_emission.dat',
                     'dust_scat.dat',
                     ext_tags=['band','iue','fuse'])
        new_DG = DustGrains()
        new_DG.from_object(DG, OD)
        DG = new_DG

    # setup the plots
    fontsize = 12
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = pyplot.subplots(ncols=3, nrows=2, figsize=(15,10))

    ws_indxs = np.argsort(DG.wavelengths)
    ews_indxs = np.argsort(DG.wavelengths_emission)
    waves = DG.wavelengths[ws_indxs]
    for i in range(DG.n_sizes):
        pcolor = colorsys.hsv_to_rgb(float(i) / DG.n_sizes / (1.1), 1, 1)

        ax[0,0].plot(waves, DG.cabs[i,ws_indxs],color=pcolor)
        ax[0,0].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[0,0].set_ylabel('C(abs)')
        ax[0,0].set_xscale('log')
        ax[0,0].set_yscale('log')

        ax[0,1].plot(waves, DG.csca[i,ws_indxs],color=pcolor)
        ax[0,1].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[0,1].set_ylabel('C(sca)')
        ax[0,1].set_xscale('log')
        ax[0,1].set_yscale('log')

        ax[0,2].plot(waves, DG.cext[i,ws_indxs],color=pcolor)
        ax[0,2].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[0,2].set_ylabel('C(sca)')
        ax[0,2].set_xscale('log')
        ax[0,2].set_yscale('log')

        ax[1,0].plot(DG.wavelengths_scat_a, 
                     DG.scat_a_csca[i,:]/DG.scat_a_cext[i,:],
                     'o',
                     color=pcolor)
        ax[1,0].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[1,0].set_ylabel('albedo')
        ax[1,0].set_xscale('log')

        ax[1,1].plot(DG.wavelengths_scat_g, DG.scat_g[i,:],'o',color=pcolor)
        ax[1,1].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[1,1].set_ylabel('g')
        ax[1,1].set_xscale('log')

        ax[1,2].plot(DG.wavelengths_emission[ews_indxs], 
                     DG.emission[i,ews_indxs],color=pcolor)
        ax[1,2].set_xlabel(r'$\lambda$ [$\mu m$]')
        ax[1,2].set_ylabel('Emission')
        ax[1,2].set_xscale('log')
        ax[1,2].set_yscale('log')
        cur_ylim = ax[1,2].get_ylim()
        ax[1,2].set_ylim([1e-23,1e-0])
        
    pyplot.tight_layout()    

    # show or save
    basename = 'DustGrains_diag'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        pyplot.show()
    
    
