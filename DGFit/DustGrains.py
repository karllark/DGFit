#!/usr/bin/env python
#Started: Jan 2015 (KDG)
#Updated to include better diagnoistic plots when run (Mar 2016 KDG) 
"""
DustGrains class
  dust grain properties stored by dust size/composition
"""
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

import matplotlib.pyplot as plt
import matplotlib

from scipy.interpolate import interp1d

from ObsData import ObsData

# Object for the proprerties of dust grain with a specific composition
class DustGrains():
    """
    DustGrains Class

    Parameters
    ----------
    
    Attributes
    ----------
    origin : 'string'

    """
    def __init__(self):
        """
        Simple initialization allowing for multiple origins of data
        """
        self.origin = None

    def from_files(self, componentname, path='./'):
        """
        Read in precomputed dust grain information from files.
        
        Parameters
        ----------
        componentname : 'string'
            Name that givesn the dust composition
            [astro-silicates, astro-carbonacenous, astro-graphite]

        path : 'string'
            Path to the location of the dust grain files
        """
        
        self.origin = 'files'
        
        # min/max wavelengths for storage
        #    set here in case later we want to pass them via the function call
        min_wave = 0.
        max_wave = 1e6,
        min_wave_emission = 0.
        max_wave_emission = 1e6

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

        # check if any files were found
        if len(filelist) == 0:
            print("no files found")
            print("path = " + path)
            exit()

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
                           header_start=-1)
            
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
            
            # default size distributions
            self.size_dist[k] = self.sizes[k]**(-4.0)

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


    def from_object(self, DustGrain, ObsData):
        """
        Setup a new DustGrains object on the ObsData object wavelength grids
        using an existing DustGrain object for the dust grain information.
        Currently the information is interpolated to the new wavelength grids.

        In the future, this should be enhanced to integrate across filter
        bandpasses for the data derived in filters.
        
        Parameters
        ----------
        DustGrain : DustGrains object 
           usually read from the files with the from_files function

        ObsData: ObsData object
           contains all the observed data to be fit
        """
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
        self.n_wavelengths = len(self.wavelengths)

        # variables to store the dust grain properties
        self.cext = np.empty((self.n_sizes, self.n_wavelengths))
        self.cabs = np.empty((self.n_sizes, self.n_wavelengths))
        self.csca = np.empty((self.n_sizes, self.n_wavelengths))

        if ObsData.fit_ir_emission:
            self.wavelengths_emission = ObsData.ir_emission_waves
            self.n_wavelengths_emission = len(self.wavelengths_emission)
            self.emission = np.empty((self.n_sizes,
                                      self.n_wavelengths_emission))

        if ObsData.fit_scat_a:
            self.wavelengths_scat_a = ObsData.scat_a_waves
            self.n_wavelengths_scat_a = len(self.wavelengths_scat_a)
            self.scat_a_cext = np.empty((self.n_sizes,
                                         self.n_wavelengths_scat_a))
            self.scat_a_csca = np.empty((self.n_sizes,
                                         self.n_wavelengths_scat_a))

        if ObsData.fit_scat_a:
            self.wavelengths_scat_g = ObsData.scat_g_waves
            self.n_wavelengths_scat_g = len(self.wavelengths_scat_g)
            self.scat_g = np.empty((self.n_sizes, self.n_wavelengths_scat_g))
            self.scat_g_csca = np.empty((self.n_sizes,
                                         self.n_wavelengths_scat_g))

        # loop over the sizes and generate grain info on the observed data grid
        for i in range(self.n_sizes):
            cext_interp = interp1d(DustGrain.wavelengths, DustGrain.cext[i,:])
            cabs_interp = interp1d(DustGrain.wavelengths, DustGrain.cabs[i,:])
            csca_interp = interp1d(DustGrain.wavelengths, DustGrain.csca[i,:])
            self.cext[i,:] = cext_interp(self.wavelengths)
            self.cabs[i,:] = cabs_interp(self.wavelengths)
            self.csca[i,:] = csca_interp(self.wavelengths)

            if ObsData.fit_scat_a:
                self.scat_a_cext[i,:] = cext_interp(self.wavelengths_scat_a)
                self.scat_a_csca[i,:] = csca_interp(self.wavelengths_scat_a)

            if ObsData.fit_scat_g:
                g_interp = interp1d(DustGrain.wavelengths, DustGrain.scat_g[i,:])
                self.scat_g[i,:] = g_interp(self.wavelengths_scat_g)
                self.scat_g_csca[i,:] = csca_interp(self.wavelengths_scat_g)

            if ObsData.fit_ir_emission:
                emission_interp = interp1d(DustGrain.wavelengths_emission, 
                                           DustGrain.emission[i,:])
                self.emission[i,:] = emission_interp(self.wavelengths_emission)

    # function to integrate this component
    # returns the effective/total cabs, csca, etc.
    # these are normalized to NHI (assumption)
    def eff_grain_props(self, ObsData):
        """
        Calculate the grain properties integrated over the size distribution
        for a single grain composition.

        Returns
        -------
        A dictonary of:

        C(abs) : 'numpy.ndarray'
           Absorption cross section

        C(sca) : 'numpy.ndarray'
           Scattering cross section

        Abundances : ('list', 'numpy.ndarray')
           Tuple with (atomic elements, # per/10^6 H atoms

        Emission : 'numpy.ndarray'
           IR emission

        albedo : 'numpy.ndarray'
           Dust scattering albedo [Albedo C(sca)/Albedo C(ext)]

        g : 'numpy.ndarray'
           Dust scattering phase function assymetry [g = <cos theta>]

        Albedo C(ext) : 'numpy.ndarray'
           Extinction cross section on the albedo wavelength grid
           (needed for combining with other dust grain compositions)

        Albedo C(sca) : 'numpy.ndarray'
           Scattering cross section on the albedo wavelength grid
           (needed for combining with other dust grain compositions)

        G C(sca) : 'numpy.ndarray'
           Scattering cross section on the g wavelength grid 
           (needed for combining with other dust grain compositions)
        """

        # output is a dictonary
        results = {}
        
        # initialize the results
        _effcabs = np.empty(self.n_wavelengths)
        _effcsca = np.empty(self.n_wavelengths)
        
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
            
        results['cabs'] = _effcabs
        results['csca'] = _effcsca
        #results = [_effcabs, _effcsca]

        # compute the number of atoms/NHI 
        _natoms = np.empty(len(self.atomic_comp_names))
        for i in range(len(self.atomic_comp_names)):
            _natoms[i] = np.sum(deltas*( ((self.sizes[0:self.n_sizes-1]**3)*
                                          self.size_dist[0:self.n_sizes-1]*
                                          self.col_den_constant[i]) +
                                         ((self.sizes[1:self.n_sizes]**3)*
                                          self.size_dist[1:self.n_sizes]*
                                          self.col_den_constant[i]) ) )

        # convert to N(N) per 1e6 N(HI)
        _natoms *= 1e6
        
        results['natoms'] = dict(zip(self.atomic_comp_names, _natoms))
        #results.append(dict(zip(self.atomic_comp_names, _natoms)))
        
        # compute the integrated emission spectrum
        if ObsData.fit_ir_emission:
            _emission = np.empty(self.n_wavelengths_emission)
            for i in range(self.n_wavelengths_emission):
                _emission[i] = np.sum( deltas*(
                    (self.emission[0:self.n_sizes-1,i]*
                     sizedist1) +
                    (self.emission[1:self.n_sizes,i]*
                     sizedist2) ) )
            results['emission'] = _emission    
            #results.append(_emission)
            
        # scattering parameters a & g
        if ObsData.fit_scat_a:
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

            results['albedo'] = _effscat_a_csca/_effscat_a_cext
            results['scat_a_cext'] = _effscat_a_cext
            results['scat_a_csca'] = _effscat_a_csca

        if ObsData.fit_scat_g:
            n_waves_scat_g = self.n_wavelengths_scat_g
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
            results['g'] = _effg/_effscat_g_csca
            results['scat_g_csca'] = _effscat_g_csca

        # return the results as a tuple of arrays
        return results

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
                     'data_mw_rv31/MW_diffuse_Gordon09_avnhi.dat',
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

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(15,10))

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
        
    plt.tight_layout()    

    # show or save
    basename = 'DustGrains_diag'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        plt.show()
    
    
