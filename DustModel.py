#!/usr/bin/env python
#
# dustmodel_props object
#  dust model properites
#
# Started: Jan 2015 (KDG)
# 
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as pyplot
from astropy.io import fits

import DustGrains

# Object for the proprerties of dust grain with a specific composition
class DustModel():
    def __init__(self, componentnames, path='./', min_wave=0., max_wave=1e6,
                 min_wave_emission=0., max_wave_emission=1e6):
        self.n_components = len(componentnames)
        self.components = []
        # get the basic grain data
        for componentname in componentnames:
            self.components.append(DustGrains.DustGrains(componentname,path=path,min_wave=min_wave,max_wave=max_wave,
                                                         min_wave_emission=min_wave_emission, max_wave_emission=max_wave_emission))

    # set the size distributions
    # new_size_dists are the concatenated size distributions for all the components
    def set_size_dist(self, new_size_dists):
        k1 = 0
        for component in self.components:
            k2 = k1 + component.n_sizes
            #print(k1,k2, component.n_sizes)
            component.size_dist[:] = new_size_dists[k1:k2]
            k1 += component.n_sizes
                
    # compute integrated dust properties
    def eff_grain_props(self):
        _cabs = np.zeros(self.components[0].n_wavelengths)
        _csca = np.zeros(self.components[0].n_wavelengths)
        _emission = np.zeros(self.components[0].n_wavelengths_emission)
        _natoms = {}
        for component in self.components:
            _tcabs, _tcsca, _tnatoms, _temission = component.eff_grain_props()
            # add the component info to the total values
            _cabs += _tcabs
            _csca += _tcsca
            _emission += _temission
            # for the depletions (# of atoms), a bit more careful work needed
            for aname in _tnatoms.keys():
                #if (len(_natoms) > 0) & (aname in _natoms.keys()):
                if aname in _natoms.keys():
                    _natoms[aname] += _tnatoms[aname]
                else:
                    _natoms[aname] = _tnatoms[aname]

        return (_cabs, _csca, _natoms, _emission)

    def save(self, filename, size_dist_uncs=[0]):

        # write a small primary header
        pheader = fits.Header()
        pheader.set('NCOMPS', len(self.components), 'number of dust grain components')
        for k, component in enumerate(self.components):
            pheader.set('CNAME'+str(k), component.name, 'name of dust grain component')
        pheader.add_comment('Dust Model reuslts written by DustModel.py')
        pheader.add_comment('written by Karl D. Gordon')
        pheader.add_comment('kgordon@stsci.edu')
        phdu = fits.PrimaryHDU(header=pheader)
        
        hdulist = fits.HDUList([phdu])

        # output the dust grain size distribution
        k1 = 0
        for component in self.components:
            col1 = fits.Column(name='SIZE', format='E', array=component.sizes)
            col2 = fits.Column(name='DIST', format='E', array=component.size_dist)
            all_cols = [col1, col2]

            k2 = k1 + component.n_sizes
            if len(size_dist_uncs) > 1:
                col3 = fits.Column(name='DISTPUNC', format='E', array=size_dist_uncs[0][k1:k2])
                all_cols.append(col3)
                col4 = fits.Column(name='DISTMUNC', format='E', array=size_dist_uncs[1][k1:k2])
                all_cols.append(col4)
            k1 += component.n_sizes

            cols = fits.ColDefs(all_cols)
            tbhdu = fits.BinTableHDU.from_columns(all_cols)
            tbhdu.header.set('EXTNAME', component.name, 'dust grain component name')

            hdulist.append(tbhdu)

        # output the resulting observable parameters
        cabs, csca, natoms, emission = self.eff_grain_props()

        # natoms
        col1 = fits.Column(name='NAME', format='A2', array=np.array(list(natoms.keys())))
        col2 = fits.Column(name='ABUND', format='E', array=np.array(list(natoms.values())))
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Abundances', 'abundances in units of # atoms/1e6 H atoms')
        hdulist.append(tbhdu)

        # extinction
        col1 = fits.Column(name='WAVE', format='E', array=self.components[0].wavelengths)
        col2 = fits.Column(name='EXT', format='E', array=1.086*(cabs+csca))
        all_cols_ext = [col1, col2]
        
        # emission
        col1 = fits.Column(name='WAVE', format='E', array=self.components[0].wavelengths_emission)
        col2 = fits.Column(name='EMIS', format='E', array=emission)
        all_cols_emis = [col1, col2]

        # albedo
        col1 = fits.Column(name='WAVE', format='E', array=self.components[0].wavelengths)
        col2 = fits.Column(name='ALBEDO', format='E', array=csca/(cabs+csca))
        all_cols_albedo = [col1, col2]

        for k, component in enumerate(self.components):
            tcabs, tcsca, tnatoms, temission = component.eff_grain_props()

            tcol = fits.Column(name='EXT'+str(k+1), format='E', array=1.086*(tcabs+tcsca))
            all_cols_ext.append(tcol)

            tcol = fits.Column(name='EMIS'+str(k+1), format='E', array=temission)
            all_cols_emis.append(tcol)

            tcol = fits.Column(name='ALBEDO'+str(k+1), format='E', array=tcsca/(tcsca+tcabs))
            all_cols_albedo.append(tcol)

        # now output the results
        #    extinction
        cols = fits.ColDefs(all_cols_ext)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Extinction', 'extinction in A(lambda)/N(HI) units')
        hdulist.append(tbhdu)

        #    emission
        cols = fits.ColDefs(all_cols_emis)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Emission', 'emission MJy/sr/H atom units')
        hdulist.append(tbhdu)

        #    albedo
        cols = fits.ColDefs(all_cols_albedo)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Albedo', 'dust scattering albedo')
        hdulist.append(tbhdu)

        hdulist.writeto(filename, clobber=True)
    
if __name__ == "__main__":
    
    #dustmodel = DustModel(['astro-silicates','astro-graphite'])
    dustmodel = DustModel(['astro-silicates','astro-carbonaceous'], path='/home/kgordon/Dirty_v2/write_grain/indiv_grain/')

    cabs, csca, natoms = dustmodel.eff_grain_props()
    
    print(natoms)

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(dustmodel.components[0].wavelengths, cabs+csca, 'g-')
    ax.plot(dustmodel.components[0].wavelengths, cabs, 'r-')
    ax.plot(dustmodel.components[0].wavelengths, csca, 'b-')
    ax.plot(dustmodel.components[0].wavelengths, csca/(cabs+csca), 'c-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    pyplot.show()
    exit()

