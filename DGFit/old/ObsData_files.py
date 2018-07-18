#!/usr/bin/env python
#
# DGFit Observed Data Object
#
# Started: Jan 2015 (KDG)
#    Revised to read in all the data from files: Mar 2016 (KDG)
#       observed data defines the wavelength grids to fit on
# 
from __future__ import print_function

import string

import numpy as np
import matplotlib.pyplot as pyplot
from astropy.table import Table

from Ext import f99

# Object for the observed dust data
class ObsData():

    # read in the data from files
    def __init__(self, ext_filenames, dep_filename, ir_emis_filename,
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
                                            
        # normalization from A(V) to N(HI)
        #   need to include this as input so that it can change
        self.Rv = 3.1
        # number is N(HI)/E(B-V) from Bohlin, Savage, & Drake (1978)
        self.avnhi = self.Rv/5.8e21   
        self.avnhi_unc = 0.0   # update once the true unc is known
        
        # change the extinction normalization from A(V) to N(HI)
        self.ext_alnhi = self.ext_alav*self.avnhi
        self.ext_alnhi_unc = np.square(self.ext_alav_unc/self.ext_alav) + \
                             np.square(self.avnhi_unc/self.avnhi)
        self.ext_alnhi_unc = self.ext_alnhi*np.sqrt(self.ext_alnhi_unc)
        
        # dust abundances
        t = Table.read(dep_filename,format='ascii.commented_header')
        self.abundance = {}
        self.total_abundance = {}
        for i in range(len(t)):
            self.abundance[t['atom'][i]] = (t['abund'][i], t['abund_unc'][i])
            self.total_abundance[t['atom'][i]] = (t['total_abund'][i],
                                                  t['total_abund_unc'][i])

        # diffuse IR emission spectrum (Gordon et al. 2014)
        t = Table.read(ir_emis_filename,format='ascii.commented_header')
        self.ir_emission_waves = np.array(t['wave'])
        self.ir_emission = np.array(t['emission'])
        self.ir_emission_uncs = np.array(t['emission_unc'])
                                         
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

        self.scat_waves = np.array(scat_waves)*1e-4
        self.scat_albedo = np.array(scat_albedo)
        self.scat_albedo_unc = np.array(scat_albedo_unc)
        self.scat_g = np.array(scat_g)
        self.scat_g_unc = np.array(scat_g_unc)
        self.scat_ref = scat_ref

        self.fit_depletions = True
        self.fit_ir_emission = True
        self.fit_scat_param = True
        
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

    testobs = ObsData(['data_mw_rv31/MW_diffuse_Gordon09_band_ext.dat',
                       'data_mw_rv31/MW_diffuse_Gordon09_iue_ext.dat',
                       'data_mw_rv31/MW_diffuse_Gordon09_fuse_ext.dat'],
                      'data_mw_rv31/MW_diffuse_Jenkins09_abundances.dat',
                      'data_mw_rv31/MW_diffuse_Gordon14_ir_emission.dat',
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

    for i in range(DG.n_sizes):
        ax[0,0].plot(DG.wavelengths, DG.cabs[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[0,0].set_xlabel(r'wavelength [$\mu$m]')
        ax[0,0].set_ylabel('C(abs)')
        ax[0,0].set_xscale('log')
        ax[0,0].set_yscale('log')

        ax[0,1].plot(DG.wavelengths, DG.csca[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[0,1].set_xlabel(r'wavelength [$\mu$m]')
        ax[0,1].set_ylabel('C(sca)')
        ax[0,1].set_xscale('log')
        ax[0,1].set_yscale('log')

        ax[0,2].plot(DG.wavelengths, DG.cext[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[0,2].set_xlabel(r'wavelength [$\mu$m]')
        ax[0,2].set_ylabel('C(sca)')
        ax[0,2].set_xscale('log')
        ax[0,2].set_yscale('log')

        ax[1,0].plot(DG.wavelengths, DG.csca[i,:]/DG.cext[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[1,0].set_xlabel(r'wavelength [$\mu$m]')
        ax[1,0].set_ylabel('albedo')
        ax[1,0].set_xscale('log')

        ax[1,1].plot(DG.wavelengths, DG.g[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[1,1].set_xlabel(r'wavelength [$\mu$m]')
        ax[1,1].set_ylabel('g')
        ax[1,1].set_xscale('log')

        ax[1,2].plot(DG.wavelengths, DG.emission[i,:],
                     color=colorsys.hsv_to_rgb(i / DG.n_sizes / (1.1), 1, 1))
        ax[1,2].set_xlabel(r'wavelength [$\mu$m]')
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
    
    
