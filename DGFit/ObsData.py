#!/usr/bin/env python
# Started: Jan 2015 (KDG)
#    Revised to read in all the data from files: Mar 2016 (KDG)
#       observed data defines the wavelength grids to fit on
"""
ObsData class
  observed data that will be used to constrain the dust model
"""
from __future__ import print_function

import numpy as np
from astropy.table import Table
from astropy.io import fits

__all__ = ["ObsData"]


# Object for the observed dust data
class ObsData():
    """
    ObsData Class

    Parameters
    ----------
    ext_filenames: list of 'string'
        filenames with the observed extincction curve

    avnhi_filenames: list of 'string'
        filename with the observed A(V)/N(HI) value + unc

    abund_filename: 'string'
        filename with the observed atomic abundances

    ir_emis_filename: 'string'
        filename with the observed infrared dust emission

    dust_scat_filename: 'string'
        filename with the observed dust scattering (a, g) parameters
        [currently not used - hard coded for MW diffuse - need to change]

    ext_tags : list of 'string'
        list of tags identifying the origin of the
        dust extinction curve segments

    Attributes
    ----------
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
    def __init__(self, ext_filenames, avnhi_filename,
                 abund_filename, ir_emis_filename,
                 dust_scat_filename, ext_tags=None):

        # extinction curve
        self.fit_extinction = True
        self.ext_waves = np.empty((0))
        self.ext_alav = np.empty((0))
        self.ext_alav_unc = np.empty((0))
        self.ext_tags = []

        if isinstance(ext_filenames, (list, tuple)):
            for i, filename in enumerate(ext_filenames):
                t = Table.read(filename, format='ascii.commented_header')
                self.ext_waves = np.concatenate([self.ext_waves,
                                                 1.0/t['wave']])
                self.ext_alav = np.concatenate([self.ext_alav,
                                                t['A(l)/A(V)']])
                self.ext_alav_unc = np.concatenate([self.ext_alav_unc,
                                                    t['unc']])
                if ext_tags is not None:
                    cur_tag = ext_tags[i]
                else:
                    cur_tag = 'Tag' + str(i+1)
                self.ext_tags = self.ext_tags + len(t['wave'])*[cur_tag]
        else:
            # assume it is a FITS file (need to add checks)
            hdulist = fits.open(ext_filenames)
            for i in range(1, len(hdulist)):
                t = hdulist[i].data
                # hack to get AzV 215 to work
                #  need to get a better file format for FITS extinction curves
                #  units, etc.
                trv = 3.65
                ext = (t['EXT']/trv) + 1
                ext_unc = t['UNC']/trv

                # only keep positive measurements
                gindxs, = np.where(ext > 0.0)
                self.ext_waves = np.concatenate([self.ext_waves,
                                                 t['WAVELENGTH'][gindxs]])
                self.ext_alav = np.concatenate([self.ext_alav, ext[gindxs]])
                self.ext_alav_unc = np.concatenate([self.ext_alav_unc,
                                                    ext_unc[gindxs]])
                self.ext_tags = self.ext_tags + \
                    len(t['WAVELENGTH'])*[hdulist[i].header['EXTNAME']]

            hdulist.close()

        # sort
        sindxs = np.argsort(self.ext_waves)
        self.ext_waves = self.ext_waves[sindxs]
        self.ext_alav = self.ext_alav[sindxs]
        self.ext_alav_unc = self.ext_alav_unc[sindxs]
        self.ext_tags = np.array(self.ext_tags)[sindxs]

        # normalization from A(V) to N(HI)
        t = Table.read(avnhi_filename,
                       format='ascii.commented_header',
                       header_start=-1)
        self.avnhi = t['Av_to_NHI'][0]
        self.avnhi_unc = t['unc'][0]

        # change the extinction normalization from A(V) to N(HI)
        self.ext_alnhi = self.ext_alav*self.avnhi
        self.ext_alnhi_unc = (np.square(self.ext_alav_unc/self.ext_alav)
                              + np.square(self.avnhi_unc/self.avnhi))
        self.ext_alnhi_unc = self.ext_alnhi*np.sqrt(self.ext_alnhi_unc)

        # dust abundances
        self.fit_abundance = False
        if abund_filename is not None:
            self.fit_abundance = True
            t = Table.read(abund_filename, format='ascii.commented_header')
            self.abundance = {}
            self.total_abundance = {}
            for i in range(len(t)):
                self.abundance[t['atom'][i]] = (t['abund'][i],
                                                t['abund_unc'][i])
                self.total_abundance[t['atom'][i]] = (t['total_abund'][i],
                                                      t['total_abund_unc'][i])

        # diffuse IR emission spectrum
        self.fit_ir_emission = False
        if ir_emis_filename is not None:
            self.fit_ir_emission = True
            t = Table.read(ir_emis_filename, format='ascii.commented_header')
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
        self.fit_scat_a = False
        self.fit_scat_g = False
        if dust_scat_filename is not None:
            self.fit_scat_a = True
            self.fit_scat_g = True
            files_dgl = ["mathis73", "morgan76", "lillie76", "toller81",
                         "murthy93", "murthy95", "petersohn97", "witt97",
                         "schiminovich01", "shalima04", "sujatha05",
                         "sujatha07", "sujatha10"]
            scat_path = "data/mw_rv31/Scat_Data/"

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

                t = Table.read(scat_path+sfile+'.dat',
                               format='ascii',
                               header_start=1)
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
