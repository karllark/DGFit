#!/usr/bin/env python
# Started: Jan 2015 (KDG)
#    Revised to read in all the data from files: Mar 2016 (KDG)
#       observed data defines the wavelength grids to fit on
"""
ObsData class
  observed data that will be used to constrain the dust model
"""
from __future__ import print_function

import string
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from astropy.table import Table

# Object for the observed dust data
class ObsData():
    """
    ObsData Class

    Parameters
    ----------
    ext_filenames: list of 'string'
        filenames with the observed extincction curve

    abund_filename: 'string'
        filename with the observed atomic abundances

    ir_emis_filename: 'string'
        filename with the observed infrared dust emission

    dust_scat_filename: 'string'
        filename with the observed dust scattering (a, g) parameters
        [currently not used - hard coded for MW diffuse - need to change]

    ext_tags : list of 'string'
        list of tags identifying the origin of the dust extinction curve segments

    Attributes
    ----------
    Rv : 'float'
        R(V) = A(V)/E(B-V) for extinction curve
    
    alnhi : float
        A(lamda)/N(HI) value for extinction curve

    alnhi_unc : float
        uncertainty in A(lamda)/N(HI) value for extinction curve

    ext_waves : 'numpy.ndarray'
        wavelengths for the extinction curve
    
    ext_alav : 'numpy.ndarray'
        extinction curve in A(lambda)/A(V) units
    
    ext_alav_unc : 'numpy.ndarray'
        extinction curve uncertainties in A(lambda)/A(V) units
    
    ext_alnhi : 'numpy.ndarray'
        extinction curve in A(lambda)/N(HI) units

    ext_alnhi_unc : 'numpy.ndarray'
        extinction curve uncertainties in A(lambda)/N(HI) units
    
    ext_tags : 'numpy.ndarray'
        string tags identifying the origin of the extinction curve measurement
    
    
    """
    
    # read in the data from files
    def __init__(self, ext_filenames, abund_filename, ir_emis_filename,
                 dust_scat_filename, ext_tags=None):

        # extinction curve
        self.ext_waves = np.empty((0))
        self.ext_alav = np.empty((0))
        self.ext_alav_unc = np.empty((0))
        self.ext_tags = []
        for i, filename in enumerate(ext_filenames):
            t = Table.read(filename, format='ascii.commented_header')
            self.ext_waves = np.concatenate([self.ext_waves, 1.0/t['wave']])
            self.ext_alav = np.concatenate([self.ext_alav, t['A(l)/A(V)']])
            self.ext_alav_unc = np.concatenate([self.ext_alav_unc, t['unc']])
            if ext_tags is not None:
                cur_tag = ext_tags[i]
            else:
                cur_tag = 'Tag' + str(i+1)
            self.ext_tags = self.ext_tags + len(t['wave'])*[cur_tag]

        # sort 
        sindxs = np.argsort(self.ext_waves)
        self.ext_waves = self.ext_waves[sindxs]
        self.ext_alav = self.ext_alav[sindxs]
        self.ext_alav_unc = self.ext_alav_unc[sindxs]
        self.ext_tags = np.array(self.ext_tags)[sindxs]
                                            
        # normalization from A(V) to N(HI)
        #   need to include this as input so that it can change
        #   current values from average of 40 curves
        #       with R(V)~3.1 in Gordon et al. (2009)
        self.Rv = 3.1
        self.avnhi = 5.7e-22
        self.avnhi_unc = 0.2e-22
        
        # change the extinction normalization from A(V) to N(HI)
        self.ext_alnhi = self.ext_alav*self.avnhi
        self.ext_alnhi_unc = np.square(self.ext_alav_unc/self.ext_alav) + \
                             np.square(self.avnhi_unc/self.avnhi)
        self.ext_alnhi_unc = self.ext_alnhi*np.sqrt(self.ext_alnhi_unc)
        
        # dust abundances
        t = Table.read(abund_filename,format='ascii.commented_header')
        self.abundance = {}
        self.total_abundance = {}
        for i in range(len(t)):
            self.abundance[t['atom'][i]] = (t['abund'][i], t['abund_unc'][i])
            self.total_abundance[t['atom'][i]] = (t['total_abund'][i],
                                                  t['total_abund_unc'][i])

        # diffuse IR emission spectrum (Gordon et al. 2014)
        #t = Table.read(ir_emis_filename,format='ascii.commented_header')
        #self.ir_emission_waves = np.array(t['wave'])
        #self.ir_emission = np.array(t['emission'])
        #self.ir_emission_unc = np.array(t['emission_unc'])

        #diffuse IR emission spectrum
        t = Table.read(ir_emis_filename,format='ascii.commented_header')
        self.ir_emission_waves = np.array(t['WAVE'])
        self.ir_emission = np.array(t['SPEC'])/1e20
        self.ir_emission_unc = np.array(t['ERROR'])/1e20
        # check if any uncs are zero
        gindxs, = np.where(self.ir_emission_unc == 0.0)
        if len(gindxs) > 0:
            self.ir_emission_unc[gindxs] = 0.1*self.ir_emission[gindxs]

        # sort 
        sindxs = np.argsort(self.ir_emission_waves)
        self.ir_emission_waves = self.ir_emission_waves[sindxs]
        self.ir_emission = self.ir_emission[sindxs]
        self.ir_emission_unc = self.ir_emission_unc[sindxs]

        # dust albedo (Gordon et al. AoD proceedings)
        files_dgl = ["mathis73","morgan76","lillie76","toller81", 
                     "murthy93","murthy95","petersohn97","witt97", 
                     "schiminovich01","shalima04","sujatha05","sujatha07",
                     "sujatha10"]
        scat_path = "/home/kgordon/Pro/Dust/Scat_Data/"

        scat_waves = []
        scat_albedo = []
        scat_albedo_unc = []
        scat_g = []
        scat_g_unc = []
        scat_ref = []
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

        # remove all the measurements with zero uncertainty
        gindxs, = np.where(np.array(scat_albedo_unc) > 0.0)
        self.scat_a_waves = np.array(scat_waves)[gindxs]*1e-4
        self.scat_albedo = np.array(scat_albedo)[gindxs]
        self.scat_albedo_unc = np.array(scat_albedo_unc)[gindxs]
        self.scat_a_ref = np.array(scat_ref)[gindxs]

        # sort 
        sindxs = np.argsort(self.scat_a_waves)
        self.scat_a_waves = self.scat_a_waves[sindxs]
        self.scat_albedo = self.scat_albedo[sindxs]
        self.scat_albedo_unc = self.scat_albedo_unc[sindxs]
        self.scat_a_ref = self.scat_a_ref[sindxs]

        # remove all the measurements with zero uncertainty
        gindxs, = np.where(np.array(scat_g_unc) > 0.0)
        self.scat_g_waves = np.array(scat_waves)[gindxs]*1e-4
        self.scat_g = np.array(scat_g)[gindxs]
        self.scat_g_unc = np.array(scat_g_unc)[gindxs]
        self.scat_g_ref = np.array(scat_ref)[gindxs]

        # sort 
        sindxs = np.argsort(self.scat_g_waves)
        self.scat_g_waves = self.scat_g_waves[sindxs]
        self.scat_g = self.scat_g[sindxs]
        self.scat_g_unc = self.scat_g_unc[sindxs]
        self.scat_g_ref = self.scat_g_ref[sindxs]

        # setup what can be fit
        self.fit_extinction = True
        self.fit_abundance = True
        self.fit_ir_emission = True
        self.fit_scat_a = True
        self.fit_scat_g = True
        
if __name__ == "__main__":
    
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    OD = ObsData(['data_mw_rv31/MW_diffuse_Gordon09_band_ext.dat',
                  'data_mw_rv31/MW_diffuse_Gordon09_iue_ext.dat',
                  'data_mw_rv31/MW_diffuse_Gordon09_fuse_ext.dat'],
                 'data_mw_rv31/MW_diffuse_Jenkins09_abundances.dat',
                 'data_mw_rv31/MW_diffuse_Compiegne11_ir_emission.dat',
                 'dust_scat.dat',
                 ext_tags=['band','iue','fuse'])

    # setup the plots
    fontsize = 12
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = pyplot.subplots(ncols=3, nrows=2, figsize=(15,10))

    ax[0,0].errorbar(OD.ext_waves, OD.ext_alnhi, yerr=OD.ext_alnhi_unc, fmt='o',
                     label='Extinction')
    ax[0,0].set_xlabel(r'$\lambda [\mu m]$')
    ax[0,0].set_ylabel(r'$A(\lambda)/N(HI)$')
    ax[0,0].set_xscale('log')
    ax[0,0].set_xlim(0.085,3.0)
    ax[0,0].legend()

    n_atoms = len(OD.abundance)
    aindxs = np.arange(n_atoms)
    width = 0.5
    atomnames = sorted(list(OD.abundance.keys()))

    ax[1,0].bar(aindxs+0.25*width, 
                [OD.total_abundance[x][0] for x in atomnames], 
                width, color='g', alpha=0.25,label='gas+dust')

    ax[1,0].errorbar(aindxs+0.75*width, 
                [OD.abundance[x][0] for x in atomnames],
                yerr=[OD.abundance[x][1] for x in atomnames],
                fmt='o', label='dust')

    ax[1,0].set_ylabel(r'$N(X)/[10^6 N(HI)]$', fontsize=fontsize)
    ax[1,0].set_xticks(aindxs+(0.75*width))
    ax[1,0].set_xticklabels(atomnames)
    ax[1,0].legend(loc=2)

    ax[0,1].errorbar(OD.ir_emission_waves, OD.ir_emission,
                     yerr=OD.ir_emission_unc, fmt='o',label="Emission")
    ax[0,1].set_xlabel(r'$\lambda [\mu m]$')
    ax[0,1].set_ylabel(r'$S$ $[MJy$ $sr^{-1}$ $N(HI)^{-1}]$')
    ax[0,1].set_xscale('log')
    ax[0,1].set_xlim(1.0,1.5e4)
    ax[0,1].set_yscale('log')
    ax[0,1].legend(loc=2)

    t = Table.read('data_mw_rv31/MW_diffuse_Mathis83_ISRF.dat',
                   format='ascii.commented_header')
    ax[1,1].plot(t['wave'],t['ISRF'],'-', label="ISRF")
    ax[1,1].set_xlabel(r'$\lambda [\mu m]$')
    ax[1,1].set_ylabel(r'ISRF [$ergs$ $cm^{-3}$ $s^{-1}$ $sr^{-1}$]')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    ax[1,1].set_xlim(0.09,1e1)
    ax[1,1].set_ylim(1e-2,1e2)
    ax[1,1].legend()

    ax[0,2].errorbar(OD.scat_a_waves, OD.scat_albedo,
                     yerr=OD.scat_albedo_unc, fmt='o', label='albedo')
    ax[0,2].set_xlabel(r'$\lambda [\mu m]$')
    ax[0,2].set_ylabel(r'$a$')
    ax[0,2].set_xscale('log')
    ax[0,2].set_xlim(0.085,3.0)
    ax[0,2].set_ylim(0.0,1.0)
    ax[0,2].legend()

    ax[1,2].errorbar(OD.scat_g_waves, OD.scat_g,
                     yerr=OD.scat_g_unc, fmt='o', 
                     label=r'$g = < \mathrm{cos} (\theta) >$')
    ax[1,2].set_xlabel(r'$\lambda [\mu m]$')
    ax[1,2].set_ylabel(r'$g$')
    ax[1,2].set_xscale('log')
    ax[1,2].set_xlim(0.085,3.0)
    ax[1,2].set_ylim(0.0,1.0)
    ax[1,2].legend()

    pyplot.tight_layout()    

    # show or save
    basename = 'ObsData_MW_Diffuse'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        pyplot.show()
    
    
