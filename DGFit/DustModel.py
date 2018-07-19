#!/usr/bin/env python
#
# dustmodel_props object
#  dust model properites
#
# Started: Jan 2015 (KDG)

from __future__ import print_function
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits

from DustGrains import DustGrains
from ObsData import ObsData

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
    def eff_grain_props(self, ObsData):
        # storage for results
        _cabs = np.zeros(self.components[0].n_wavelengths)
        _csca = np.zeros(self.components[0].n_wavelengths)
        _natoms = {}

        if ObsData.fit_ir_emission:
            _emission = np.zeros(self.components[0].n_wavelengths_emission)

        if ObsData.fit_scat_a:
            # _albedo = np.zeros(self.components[0].n_wavelengths_scat_a)
            _scat_a_cext = np.zeros(self.components[0].n_wavelengths_scat_a)
            _scat_a_csca = np.zeros(self.components[0].n_wavelengths_scat_a)

        if ObsData.fit_scat_g:
            _g = np.zeros(self.components[0].n_wavelengths_scat_g)
            _scat_g_csca = np.zeros(self.components[0].n_wavelengths_scat_g)

        # loop over components and accumulate the answer
        for component in self.components:
            results = component.eff_grain_props(ObsData)

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

            if ObsData.fit_ir_emission:
                _temission = results['emission']
                _emission += _temission

            if ObsData.fit_scat_a:
                _tscat_a_cext = results['scat_a_cext']
                _tscat_a_csca = results['scat_a_csca']
                _scat_a_cext += _tscat_a_cext
                _scat_a_csca += _tscat_a_csca

            if ObsData.fit_scat_g:
                _tg = results['g']
                _tscat_g_csca = results['scat_g_csca']
                _g += _tscat_g_csca*_tg
                _scat_g_csca += _tscat_g_csca

        results = {}
        results['cabs'] = _cabs
        results['csca'] = _csca
        results['natoms'] = _natoms

        if ObsData.fit_ir_emission:
            results['emission'] = _emission

        if ObsData.fit_scat_a:
            results['albedo'] = _scat_a_csca/_scat_a_cext

        if ObsData.fit_scat_g:
            results['g'] = _g/_scat_g_csca

        return results

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
        results = self.eff_grain_props(ObsData)
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
        if ObsData.fit_ir_emission:
            emission = results['emission']
            col1 = fits.Column(name='WAVE', format='E',
                               array=self.components[0].wavelengths_emission)
            col2 = fits.Column(name='EMIS', format='E',
                               array=emission)
            all_cols_emis = [col1, col2]

        # albedo
        if ObsData.fit_scat_a:
            albedo = results['albedo']
            tvals = self.components[0].wavelengths_scat_a
            col1 = fits.Column(name='WAVE', format='E',
                               array=tvals)
            col2 = fits.Column(name='ALBEDO', format='E',
                               array=albedo)
            all_cols_albedo = [col1, col2]

        # g
        if ObsData.fit_scat_g:
            g = results['g']
            tvals = self.components[0].wavelengths_scat_g
            col1 = fits.Column(name='WAVE', format='E',
                               array=tvals)
            col2 = fits.Column(name='G', format='E',
                               array=g)
            all_cols_g = [col1, col2]

        for k, component in enumerate(self.components):
            results = component.eff_grain_props(ObsData)
            tcabs = results['cabs']
            tcsca = results['csca']
            # tnatoms = results['natoms']

            tcol = fits.Column(name='EXT'+str(k+1), format='E',
                               array=1.086*(tcabs+tcsca))
            all_cols_ext.append(tcol)

            if ObsData.fit_ir_emission:
                temission = results['emission']
                tcol = fits.Column(name='EMIS'+str(k+1), format='E',
                                   array=temission)
                all_cols_emis.append(tcol)

            if ObsData.fit_scat_a:
                talbedo = results['albedo']
                tcol = fits.Column(name='ALBEDO'+str(k+1), format='E',
                                   array=talbedo)
                all_cols_albedo.append(tcol)

            if ObsData.fit_scat_g:
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
        if ObsData.fit_ir_emission:
            cols = fits.ColDefs(all_cols_emis)
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.header.set('EXTNAME', 'Emission',
                             'emission MJy/sr/H atom units')
            hdulist.append(tbhdu)

        #    albedo
        if ObsData.fit_scat_a:
            cols = fits.ColDefs(all_cols_albedo)
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.header.set('EXTNAME', 'Albedo', 'dust scattering albedo')
            hdulist.append(tbhdu)

        #    g
        if ObsData.fit_scat_g:
            cols = fits.ColDefs(all_cols_g)
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.header.set('EXTNAME', 'G',
                             'dust scattering phase function asymmetry')
            hdulist.append(tbhdu)

        hdulist.writeto(filename, clobber=True)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--obsdata", help="transform to observed data grids",
                        action="store_true")
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    OD = ObsData(['data/mw_rv31/MW_diffuse_Gordon09_band_ext.dat',
                  'data/mw_rv31/MW_diffuse_Gordon09_iue_ext.dat',
                  'data/mw_rv31/MW_diffuse_Gordon09_fuse_ext.dat'],
                 'data/mw_rv31/MW_diffuse_Gordon09_avnhi.dat',
                 'data/mw_rv31/MW_diffuse_Jenkins09_abundances.dat',
                 'data/mw_rv31/MW_diffuse_Compiegne11_ir_emission.dat',
                 'dust_scat.dat',
                 ext_tags=['band', 'iue', 'fuse'])

    # setup the plots
    fontsize = 12
    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    # dustmodel = DustModel(['astro-silicates','astro-graphite'])
    DM = DustModel()
    DM.predict_full_grid(['astro-silicates', 'astro-carbonaceous'],
                         path='data/indiv_grain/')

    if args.obsdata:
        DM_obs = DustModel()
        DM_obs.predict_observed_data(DM, OD)
        DM = DM_obs

    results = DM.eff_grain_props(OD)
    cabs = results['cabs']
    csca = results['csca']
    natoms = results['natoms']
    emission = results['emission']
    albedo = results['albedo']
    g = results['g']

    fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(16, 12))

    # plot the total results
    ax[1, 0].plot(DM.components[0].wavelengths, cabs+csca, 'k-')
    ax[1, 0].set_xscale('log')
    ax[1, 0].set_yscale('log')
    ax[1, 0].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[1, 0].set_ylabel(r'C(ext)')
    # ax[1,0].set_xlim(1e-2,1e0)
    # ax[1,0].set_ylim(1e3,1e5)

    ax[1, 1].plot(DM.components[0].wavelengths_emission, emission, 'k-')
    ax[1, 1].set_xscale('log')
    ax[1, 1].set_yscale('log')
    ax[1, 1].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[1, 1].set_ylabel(r'S')
    # ax[1,1].set_xlim(1e0,1e4)
    gindxs, = np.where(DM.components[0].wavelengths > 1e0)
    # ax[1,1].set_ylim(min(emission[gindxs]), max(emission[gindxs]))

    ax[1, 2].plot(DM.components[0].wavelengths_scat_a, albedo, 'k-')
    ax[1, 2].set_xscale('log')
    ax[1, 2].set_yscale('linear')
    ax[1, 2].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[1, 2].set_ylabel(r'$a$')

    ax[2, 2].plot(DM.components[0].wavelengths_scat_g, g, 'k-')
    ax[2, 2].set_xscale('log')
    ax[2, 2].set_yscale('linear')
    ax[2, 2].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[2, 2].set_ylabel(r'$g$')

    # plot the size distributions and component results
    for component in DM.components:
        ax[0, 0].plot(component.sizes, component.size_dist, '-',
                      label=component.name)
        ax[0, 1].plot(component.sizes,
                      np.power(component.sizes, 4.0)*component.size_dist, '-',
                      label=component.name)
        ax[0, 2].plot(component.sizes,
                      np.power(component.sizes, 3.0)*component.size_dist, '-',
                      label=component.name)

        cresults = component.eff_grain_props(OD)
        ax[1, 0].plot(component.wavelengths, cresults['cabs']+cresults['csca'])
        ax[1, 1].plot(component.wavelengths_emission, cresults['emission'])
        ax[1, 2].plot(component.wavelengths_scat_a, cresults['albedo'])
        ax[2, 2].plot(component.wavelengths_scat_g, cresults['g'])

    ax[0, 0].set_xscale('log')
    ax[0, 0].set_yscale('log')
    ax[0, 0].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[0, 0].set_ylabel(r'$N(a)$')

    ax[0, 1].set_xscale('log')
    ax[0, 1].set_yscale('log')
    ax[0, 1].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[0, 1].set_ylabel(r'$a^4 N(a)$')

    ax[0, 2].set_xscale('log')
    ax[0, 2].set_yscale('log')
    ax[0, 2].set_xlabel(r'$\lambda$ [$\mu m$]')
    ax[0, 2].set_ylabel(r'$a^3 N(a)$')

    ax[0, 0].legend()

    plt.tight_layout()

    # show or save
    basename = 'ObsData_MW_Diffuse'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        plt.show()
