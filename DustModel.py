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

    def save(self, filename):

        # write a small primary header
        pheader = fits.Header()
        pheader.add_comment('Dust Model reuslts written by DustModel.py')
        pheader.add_comment('programs written by Karl D. Gordon')
        pheader.add_comment('kgordon@stsci.edu')
        phdu = fits.PrimaryHDU(header=pheader)
        
        hdulist = fits.HDUList([phdu])

        for component in self.components:
            col1 = fits.Column(name='SIZE', format='E', array=component.sizes)
            col2 = fits.Column(name='DIST', format='E', array=component.size_dist)
            #col3 = fits.Column(name='DISTPUNC', format='E', array=component.size_dist_punc)
            #col4 = fits.Column(name='DISTMUNC', format='E', array=component.size_dist_munc)
            #cols = fits.ColDefs([col1, col2, col3, col4])
            cols = fits.ColDefs([col1, col2])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.header.set('EXTNAME', component.name, 'dust grain component name')
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

