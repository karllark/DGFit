#!/usr/bin/env python
#
# dustmodel_props object
#  dust model properites
#
# Started: Jan 2015 (KDG)

from __future__ import print_function

import numpy as np
from astropy.io import fits

from DGFit.DustGrains import DustGrains

__all__ = ["DustModel"]


# Object for the proprerties of dust grain with a specific composition
class DustModel():
    def __init__(self):
        self.origin = None

    def predict_full_grid(self, componentnames, path='./',
                          min_wave=0., max_wave=1e6,
                          min_wave_emission=0., max_wave_emission=1e6):
        self.origin = 'files'
        self.n_components = len(componentnames)
        self.components = []
        # get the basic grain data
        for componentname in componentnames:
            cur_DG = DustGrains()
            cur_DG.from_files(componentname,
                              path=path)
            self.components.append(cur_DG)

    # calculate the dust grain properties in the observed data space
    #   basically, transform the unifrom dust grain grid to the
    #   the nonuniform spectroscipic and band integrated grids
    #   of the observed data
    # this is caching the dust grains predictions to make the fitting faster
    def predict_observed_data(self, DustModel, ObsData):

        self.origin = 'obsdata'
        self.n_components = DustModel.n_components
        self.components = []
        for component in DustModel.components:
            cur_DG = DustGrains()
            cur_DG.from_object(component, ObsData)
            self.components.append(cur_DG)

    # set the size distributions
    # new_size_dists are the concatenated size distributions for
    #     all the components
    # input as a log to ensure it is always positive
    def set_size_dist(self, new_size_dists):
        k1 = 0
        for component in self.components:
            k2 = k1 + component.n_sizes
            component.size_dist[:] = new_size_dists[k1:k2]
            k1 += component.n_sizes

    # compute integrated dust properties
    def eff_grain_props(self, ObsData,
                        predict_all=False):
        # storage for results
        _cabs = np.zeros(self.components[0].n_wavelengths)
        _csca = np.zeros(self.components[0].n_wavelengths)
        _natoms = {}

        if ObsData.fit_ir_emission or predict_all:
            _emission = np.zeros(self.components[0].n_wavelengths_emission)

        if ObsData.fit_scat_a or predict_all:
            # _albedo = np.zeros(self.components[0].n_wavelengths_scat_a)
            _scat_a_cext = np.zeros(self.components[0].n_wavelengths_scat_a)
            _scat_a_csca = np.zeros(self.components[0].n_wavelengths_scat_a)

        if ObsData.fit_scat_g or predict_all:
            _g = np.zeros(self.components[0].n_wavelengths_scat_g)
            _scat_g_csca = np.zeros(self.components[0].n_wavelengths_scat_g)

        # loop over components and accumulate the answer
        for component in self.components:
            results = component.eff_grain_props(ObsData,
                                                predict_all=predict_all)

            # add the component info to the total values
            _tcabs = results['cabs']
            _tcsca = results['csca']
            _cabs += _tcabs
            _csca += _tcsca

            # for the depletions (# of atoms), a bit more careful work needed
            _tnatoms = results['natoms']
            for aname in _tnatoms.keys():
                # if (len(_natoms) > 0) & (aname in _natoms.keys()):
                if aname in _natoms.keys():
                    _natoms[aname] += _tnatoms[aname]
                else:
                    _natoms[aname] = _tnatoms[aname]

            if ObsData.fit_ir_emission or predict_all:
                _temission = results['emission']
                _emission += _temission

            if ObsData.fit_scat_a or predict_all:
                _tscat_a_cext = results['scat_a_cext']
                _tscat_a_csca = results['scat_a_csca']
                _scat_a_cext += _tscat_a_cext
                _scat_a_csca += _tscat_a_csca

            if ObsData.fit_scat_g or predict_all:
                _tg = results['g']
                _tscat_g_csca = results['scat_g_csca']
                _g += _tscat_g_csca*_tg
                _scat_g_csca += _tscat_g_csca

        results = {}
        results['cabs'] = _cabs
        results['csca'] = _csca
        results['natoms'] = _natoms

        if ObsData.fit_ir_emission or predict_all:
            results['emission'] = _emission

        if ObsData.fit_scat_a or predict_all:
            results['albedo'] = _scat_a_csca/_scat_a_cext

        if ObsData.fit_scat_g or predict_all:
            results['g'] = _g/_scat_g_csca

        return results

    def sizedist_from_file(self, filename):
        """
        Read in the size distribution from a file interpolating
        if needed

        Parameters
        ----------
        filename : str
            name of FITS file with size distributions
            one component per extension
        """
        for k, component in enumerate(self.components):
            fitsdata = fits.getdata(filename, k+1)

            # interpolate, otherwise assume exact match in sizes
            #   might want to add some checking here for robustness
            if len(component.size_dist) != len(fitsdata[:][1]):
                component.size_dist = 10**np.interp(np.log10(component.sizes),
                                                    np.log10(fitsdata['SIZE']),
                                                    np.log10(fitsdata['DIST']))
            else:
                component.size_dist = fitsdata['DIST']

    def save(self, filename, ObsData, size_dist_uncs=[0]):

        # write a small primary header
        pheader = fits.Header()
        pheader.set('NCOMPS', len(self.components),
                    'number of dust grain components')
        for k, component in enumerate(self.components):
            pheader.set('CNAME'+str(k), component.name,
                        'name of dust grain component')
        pheader.add_comment('Dust Model reuslts written by DustModel.py')
        pheader.add_comment('written by Karl D. Gordon')
        pheader.add_comment('kgordon@stsci.edu')
        phdu = fits.PrimaryHDU(header=pheader)

        hdulist = fits.HDUList([phdu])

        # output the dust grain size distribution
        k1 = 0
        for component in self.components:
            col1 = fits.Column(name='SIZE', format='E',
                               array=component.sizes)
            col2 = fits.Column(name='DIST', format='E',
                               array=component.size_dist)
            all_cols = [col1, col2]

            k2 = k1 + component.n_sizes
            if len(size_dist_uncs) > 1:
                col3 = fits.Column(name='DISTPUNC', format='E',
                                   array=size_dist_uncs[0][k1:k2])
                all_cols.append(col3)
                col4 = fits.Column(name='DISTMUNC', format='E',
                                   array=size_dist_uncs[1][k1:k2])
                all_cols.append(col4)
            k1 += component.n_sizes

            cols = fits.ColDefs(all_cols)
            tbhdu = fits.BinTableHDU.from_columns(all_cols)
            tbhdu.header.set('EXTNAME', component.name,
                             'dust grain component name')

            hdulist.append(tbhdu)

        # output the resulting observable parameters
        results = self.eff_grain_props(ObsData, predict_all=True)
        cabs = results['cabs']
        csca = results['csca']
        natoms = results['natoms']

        # natoms
        col1 = fits.Column(name='NAME', format='A2',
                           array=np.array(list(natoms.keys())))
        col2 = fits.Column(name='ABUND', format='E',
                           array=np.array(list(natoms.values())))
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Abundances',
                         'abundances in units of # atoms/1e6 H atoms')
        hdulist.append(tbhdu)

        # extinction
        col1 = fits.Column(name='WAVE', format='E',
                           array=self.components[0].wavelengths)
        col2 = fits.Column(name='EXT', format='E',
                           array=1.086*(cabs+csca))
        all_cols_ext = [col1, col2]

        # emission
        # if ObsData.fit_ir_emission:
        emission = results['emission']
        col1 = fits.Column(name='WAVE', format='E',
                           array=self.components[0].wavelengths_emission)
        col2 = fits.Column(name='EMIS', format='E',
                           array=emission)
        all_cols_emis = [col1, col2]

        # albedo
        albedo = results['albedo']
        tvals = self.components[0].wavelengths_scat_a
        col1 = fits.Column(name='WAVE', format='E',
                           array=tvals)
        col2 = fits.Column(name='ALBEDO', format='E',
                           array=albedo)
        all_cols_albedo = [col1, col2]

        # g
        g = results['g']
        tvals = self.components[0].wavelengths_scat_g
        col1 = fits.Column(name='WAVE', format='E',
                           array=tvals)
        col2 = fits.Column(name='G', format='E',
                           array=g)
        all_cols_g = [col1, col2]

        for k, component in enumerate(self.components):
            results = component.eff_grain_props(ObsData, predict_all=True)
            tcabs = results['cabs']
            tcsca = results['csca']
            # tnatoms = results['natoms']

            tcol = fits.Column(name='EXT'+str(k+1), format='E',
                               array=1.086*(tcabs+tcsca))
            all_cols_ext.append(tcol)

            temission = results['emission']
            tcol = fits.Column(name='EMIS'+str(k+1), format='E',
                               array=temission)
            all_cols_emis.append(tcol)

            talbedo = results['albedo']
            tcol = fits.Column(name='ALBEDO'+str(k+1), format='E',
                               array=talbedo)
            all_cols_albedo.append(tcol)

            tg = results['g']
            tcol = fits.Column(name='G'+str(k+1), format='E',
                               array=tg)
            all_cols_g.append(tcol)

        # now output the results
        #    extinction
        cols = fits.ColDefs(all_cols_ext)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Extinction',
                         'extinction in A(lambda)/N(HI) units')
        hdulist.append(tbhdu)

        #    emission
        cols = fits.ColDefs(all_cols_emis)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Emission',
                         'emission MJy/sr/H atom units')
        hdulist.append(tbhdu)

        #    albedo
        cols = fits.ColDefs(all_cols_albedo)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'Albedo', 'dust scattering albedo')
        hdulist.append(tbhdu)

        #    g
        cols = fits.ColDefs(all_cols_g)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.header.set('EXTNAME', 'G',
                         'dust scattering phase function asymmetry')
        hdulist.append(tbhdu)

        hdulist.writeto(filename, overwrite=True)
