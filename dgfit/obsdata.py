import numpy as np
from astropy.table import Table

__all__ = ["ObsData"]


# Object for the observed dust data
class ObsData(object):
    """
    ObsData Class

    observed data that will be used to constrain the dust model

    Parameters
    ----------
    obs_filename: 'string'
        filename with the observed data types and filenames

    Attributes
    ----------
    alnhi : float
        A(lambda)/N(HI) value for extinction curve

    alnhi_unc : float
        uncertainty in A(lambda)/N(HI) value for extinction curve

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

    """

    # read in the data from files
    def __init__(
        self,
        obs_filename,
        path = './'
    ):

        if path != './':
            obs_filename = path + obs_filename

        # get the observed data filenames
        self.parse_obsfile(obs_filename)

        # extinction curve
        self.fit_extinction = False
        if self.obs_filenames["ext"] is not None:
            self.fit_extinction = True
            t = Table.read(path + self.obs_filenames["ext"], format="ascii.commented_header")
            self.ext_waves = np.array(t["wave"])
            self.ext_alav = np.array(t["A(l)/A(V)"])
            self.ext_alav_unc = np.array(t["unc"])
            if "type" in t.colnames:
                self.ext_type = np.array(t["type"])
            else:
                self.ext_type = np.full(len(t), "spec")

            # sort
            sindxs = np.argsort(self.ext_waves)
            self.ext_waves = self.ext_waves[sindxs]
            self.ext_alav = self.ext_alav[sindxs]
            self.ext_alav_unc = self.ext_alav_unc[sindxs]
            self.ext_type = self.ext_type[sindxs]
        else:
            self.ext_waves = np.logspace(np.log10(0.0912), np.log10(32.0), 200)

        # normalization from A(V) to N(HI)
        if self.obs_filenames["avnhi"] is not None:
            t = Table.read(
                path + self.obs_filenames["avnhi"],
                format="ascii.commented_header",
                header_start=-1,
            )
            self.avnhi = t["Av_to_NHI"][0]
            self.avnhi_unc = t["unc"][0]

            # change the extinction normalization from A(V) to N(HI)
            self.ext_alnhi = self.ext_alav * self.avnhi
            self.ext_alnhi_unc = np.square(
                self.ext_alav_unc / self.ext_alav
            ) + np.square(self.avnhi_unc / self.avnhi)
            self.ext_alnhi_unc = self.ext_alnhi * np.sqrt(self.ext_alnhi_unc)
            self.ext_alnhi_npts = len(self.ext_alnhi)

        # dust abundances
        self.fit_abundance = False
        if self.obs_filenames["abund"] is not None:
            self.fit_abundance = True
            t = Table.read(path + self.obs_filenames["abund"], format="ascii.commented_header")
            self.abundance = {}
            self.total_abundance = {}
            for i in range(len(t)):
                self.abundance[t["atom"][i]] = (t["abund"][i], t["abund_unc"][i])
                self.total_abundance[t["atom"][i]] = (
                    t["total_abund"][i],
                    t["total_abund_unc"][i],
                )
            self.abundance_npts = len(self.abundance)

        # diffuse IR emission spectrum
        self.fit_ir_emission = False
        if self.obs_filenames["ir_emis"] is not None:
            self.fit_ir_emission = True
            t = Table.read(
                path + self.obs_filenames["ir_emis"], format="ascii.commented_header"
            )
            self.ir_emission_waves = np.array(t["WAVE"])
            self.ir_emission = np.array(t["SPEC"]) / 1e20
            self.ir_emission_unc = np.array(t["ERROR"]) / 1e20
            # check if any uncs are zero
            (gindxs,) = np.where(self.ir_emission_unc == 0.0)
            if len(gindxs) > 0:
                self.ir_emission_unc[gindxs] = 0.1 * self.ir_emission[gindxs]

            # sort
            sindxs = np.argsort(self.ir_emission_waves)
            self.ir_emission_waves = self.ir_emission_waves[sindxs]
            self.ir_emission = self.ir_emission[sindxs]
            self.ir_emission_unc = self.ir_emission_unc[sindxs]
        else:
            self.ir_emission_waves = np.logspace(np.log10(1.0), np.log10(500.0), 100)

        # dust albedo (Gordon et al. 2004 AoD proceedings)
        self.fit_scat_a = False
        self.fit_scat_g = False
        if (self.obs_filenames["scat_a"] is not None) & (
            self.obs_filenames["scat_g"] is not None
        ):
            self.fit_scat_a = True
            self.fit_scat_g = True

            t = Table.read(
                path + self.obs_filenames["scat_a"], format="ascii.commented_header"
            )
            self.scat_a_waves = np.array(t["wave"])
            self.scat_albedo = np.array(t["albedo"])
            self.scat_albedo_unc = np.array(t["unc"])

            t = Table.read(
                path + self.obs_filenames["scat_g"], format="ascii.commented_header"
            )
            self.scat_g_waves = np.array(t["wave"])
            self.scat_g = np.array(t["g"])
            self.scat_g_unc = np.array(t["unc"])

        else:
            self.scat_a_waves = np.logspace(np.log10(0.1), np.log10(5.0), 100)
            self.scat_g_waves = np.logspace(np.log10(0.1), np.log10(5.0), 100)

    def parse_obsfile(self, obs_filename):
        """
        Pares the observations file and get the filenames for the
        different kinds of observations
        """
        self.obs_filenames = {}

        # parse the observations file
        t = Table.read(obs_filename, format="ascii.commented_header")

        poss_types = ["ext", "ir_emis", "abund", "avnhi", "scat_a", "scat_g"]
        for ctype in poss_types:
            if ctype in t["type"]:
                gvals = t["type"] == ctype
                self.obs_filenames[ctype] = t["filename"][gvals].data[0]
            else:
                self.obs_filenames[ctype] = None
