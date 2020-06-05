import math
import numpy as np
from scipy.special import erf

from astropy.io import fits

from dgfit.dustgrains import DustGrains

__all__ = ['DustModel',
           'MRNDustModel',
           'WDDustModel']


class DustModel(object):
    """
    Full dust model including arbirary size and composition distributions.
    Includes the physical properties of the individual dust grains.

    Dust model that has each bin as an independent variable in the
    grain size distribution providing a truly arbitrary specification.

    Parameters
    ----------
    componentnames : str list, optional
        if set, then read in the grain information from files
    path : str, optional
        path to grain files
    dustmodel : DustModel object, optional
        if set, create the grain info on the obsdata wavelengths using
        the input dustmodel grain information
    obsdata : ObsData object, optional
        observed data information

    Attributes
    ----------
    origin : string
        origin of the dust grain physical properties
        allowed values are 'files' and 'onobsdata'
    n_components : int
        number of dust grain components
    components : array of DustGrain objects
        one DustGrain object per component
    sizedisttype : string
        functional form of component size distributions
    n_params : ints
        number of size distribution parameters per grain component
    parameters : dict
        Dictonary of parameters with an entry for each composition
        each entry is then a dictonary giving the value by parameter name.
        For the bins case, the dictonary is empty as the parameters is
        the size distribution.
    """
    def __init__(self,
                 componentnames=None, path='./',
                 dustmodel=None, obsdata=None,
                 every_nth=5):
        self.origin = None
        self.n_components = 0
        self.components = []
        self.sizedisttype = 'bins'
        self.n_params = None
        self.parameters = {}

        # populate the grain info
        if componentnames is not None:
            self.read_grain_files(componentnames, path=path,
                                  every_nth=every_nth)
        elif dustmodel is not None:
            self.grains_on_obs(dustmodel, obsdata)

        # set the number of size distribution parametres
        if self.n_components > 0:
            self.n_params = []
            for component in self.components:
                self.n_params.append(component.n_sizes)

    def read_grain_files(self, componentnames, path='./',
                         every_nth=5):
        """
        Read in the precomputed dust grain physical properties from files
        for each grain component.

        Parameters
        ----------
        componentnames : list of strings
            names of dust grain materials
        path : type
            path to files
        every_nth : int
            Only use every nth size, faster fitting

        Returns
        -------
        updated class variables
        """
        self.origin = 'files'
        self.n_components = len(componentnames)
        # get the basic grain data
        for componentname in componentnames:
            cur_DG = DustGrains()
            cur_DG.from_files(componentname,
                              path=path, every_nth=every_nth)
            self.components.append(cur_DG)

    def grains_on_obs(self, full_dustmodel, observeddata):
        """
        Calculate the dust grain properties on the observed
        wavelength grid.  Uses an existing DustModel based
        on the full precomputed files and an ObsData object
        to get the wavelength grid.  Makes the fitting faster
        to only do this transformation once.

        Parameters
        ----------
        full_dustmodel : DustModel object
            full dust model based on input files
        observeddata: ObsData object
            observed data to use for transformation

        Returns
        -------
        updated class variables
        """
        self.origin = 'onobsdata'
        self.n_components = full_dustmodel.n_components
        for component in full_dustmodel.components:
            cur_DG = DustGrains()
            cur_DG.from_object(component, observeddata)
            self.components.append(cur_DG)

    def compute_size_dist(self, x, params):
        """
        Compute the size distribution for the input sizes.
        For the bins case, just passes the parameters back.  Allows for
        other functional forms of the size distribution with minimal new code.

        Parameters
        ----------
        x : floats
            grains sizes
        params : floats
            Size distribution parameters
            For the arbitrary bins case, the parameters are the number
            of grains per size distribution

        Returns
        -------
        floats
            Size distribution as a function of x
        """
        return params

    def set_size_dist_parameters(self, params):
        """
        Set the size distribution parameters in the object dictonary.
        For the bins case, this does nothing.  Allows for
        other functional forms of the size distribution with minimal new code.

        Parameters
        ----------
        params : floats
            Size distribution parameters
            For the arbitrary bins case, the parameters are the number
            of grains per size distribution
        """
        pass

    def set_size_dist(self, params):
        """
        Set the size distributions for each component based on the
        parameters of the functional form of the distributions.

        Parameters
        ----------
        new_size_dists : type
            Description of parameter `new_size_dists`.

        Returns
        -------
        type
            Description of returned object.

        """
        k1 = 0
        for k, component in enumerate(self.components):
            delta_val = self.n_params[k]
            k2 = k1 + delta_val
            component.size_dist[:] = self.compute_size_dist(component.sizes[:],
                                                            params[k1:k2])
            k1 += delta_val

    def eff_grain_props(self, OD,
                        predict_all=False):
        """
        Compute the effective grain properties of the ensemble of grain
        sizes and compositions.

        Parameters
        ----------
        OD : ObsData object
            Observed data object specifically used to determine which
            observations to compute (only those needed for speed)
        predict_all : type
            Regardless of the ObsData, compute all possible observations

        Returns
        -------
        dict
            Dictonary of predicted observations
            E.g., keys of cext, natoms, emission, albedo, g
        """
        # storage for results
        _cabs = np.zeros(self.components[0].n_wavelengths)
        _csca = np.zeros(self.components[0].n_wavelengths)
        _natoms = {}

        if OD.fit_ir_emission or predict_all:
            _emission = np.zeros(self.components[0].n_wavelengths_emission)

        if OD.fit_scat_a or predict_all:
            _scat_a_cext = np.zeros(self.components[0].n_wavelengths_scat_a)
            _scat_a_csca = np.zeros(self.components[0].n_wavelengths_scat_a)

        if OD.fit_scat_g or predict_all:
            _g = np.zeros(self.components[0].n_wavelengths_scat_g)
            _scat_g_csca = np.zeros(self.components[0].n_wavelengths_scat_g)

        for component in self.components:
            results = component.eff_grain_props(OD,
                                                predict_all=predict_all)

            _tcabs = results['cabs']
            _tcsca = results['csca']
            _cabs += _tcabs
            _csca += _tcsca

            # for the depletions (# of atoms), a bit more careful work needed
            _tnatoms = results['natoms']
            for aname in _tnatoms.keys():
                if aname in _natoms.keys():
                    _natoms[aname] += _tnatoms[aname]
                else:
                    _natoms[aname] = _tnatoms[aname]

            if OD.fit_ir_emission or predict_all:
                _temission = results['emission']
                _emission += _temission

            if OD.fit_scat_a or predict_all:
                _tscat_a_cext = results['scat_a_cext']
                _tscat_a_csca = results['scat_a_csca']
                _scat_a_cext += _tscat_a_cext
                _scat_a_csca += _tscat_a_csca

            if OD.fit_scat_g or predict_all:
                _tg = results['g']
                _tscat_g_csca = results['scat_g_csca']
                _g += _tscat_g_csca*_tg
                _scat_g_csca += _tscat_g_csca

        results = {}
        results['cabs'] = _cabs
        results['csca'] = _csca
        results['natoms'] = _natoms

        if OD.fit_ir_emission or predict_all:
            results['emission'] = _emission

        if OD.fit_scat_a or predict_all:
            results['albedo'] = _scat_a_csca/_scat_a_cext

        if OD.fit_scat_g or predict_all:
            results['g'] = _g/_scat_g_csca

        return results

    def read_sizedist_from_file(self, filename):
        """
        Read in the size distribution from a file interpolating
        across sizes if needed

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

    def lnprob_generic(self, obsdata):
        """
        Compute the ln(prob) for the dust grain size and composition
        distribution as defined by the dustmodel.

        Parameters
        ----------
        obsdata : ObsData object
            All the observed data

        Returns
        -------
        float
            natural log of the probability
        """
        # get the integrated dust properties
        results = self.eff_grain_props(obsdata)

        # compute the ln(prob) for A(l)/N(HI)
        lnp_alnhi = 0.0
        if obsdata.fit_extinction:
            cabs = results['cabs']
            csca = results['csca']
            cext = cabs + csca
            dust_alnhi = 1.086*cext
            lnp_alnhi = -0.5*np.sum(((obsdata.ext_alnhi - dust_alnhi)
                                     / obsdata.ext_alnhi_unc)**2)
        # lnp_alnhi /= obsdata.n_wavelengths

        # compute the ln(prob) for the depletions
        lnp_dep = 0.0
        if obsdata.fit_abundance:
            natoms = results['natoms']
            for atomname in natoms.keys():
                lnp_dep = ((natoms[atomname] -
                            obsdata.abundance[atomname][0])
                           / obsdata.abundance[atomname][1])**2
            lnp_dep *= -0.5

        # compute the ln(prob) for IR emission
        lnp_emission = 0.0
        if obsdata.fit_ir_emission:
            emission = results['emission']
            lnp_emission = -0.5*np.sum((((obsdata.ir_emission - emission)
                                        / (obsdata.ir_emission_unc))**2))

        # compute the ln(prob) for the dust albedo
        lnp_albedo = 0.0
        if obsdata.fit_scat_a:
            albedo = results['albedo']
            lnp_albedo = -0.5*np.sum((((obsdata.scat_albedo - albedo)
                                     / (obsdata.scat_albedo_unc))**2))

        # compute the ln(prob) for the dust g
        lnp_g = 0.0
        if obsdata.fit_scat_g:
            g = results['g']
            lnp_albedo = -0.5*np.sum((((obsdata.scat_g - g)
                                       / (obsdata.scat_g_unc))**2))

        # combine the lnps
        lnp = lnp_alnhi + lnp_dep + lnp_emission + lnp_albedo + lnp_g

        # print(params)
        # print(lnp_alnhi, lnp_dep, lnp_emission, lnp_albedo, lnp_g)

        if math.isinf(lnp) | math.isnan(lnp):
            print(lnp_alnhi, lnp_dep, lnp_emission, lnp_albedo, lnp_g)
            print(lnp)
            # print(params)
            exit()
        else:
            return lnp

    @staticmethod
    def lnprob(params, obsdata, dustmodel):
        """
        Compute the full probability function including priors
        Static function as it will be called form the fitter

        Parameters
        ----------
        params : floats
            Parameters of the size distribution function
        obsdata : ObsData object
            Observed data to be fit
        dustmodel : DustModel object
            Dust model information

        Returns
        -------
        float
            natural log of the probability
        """
        # prior
        #    make sure the size distributions are all positve
        lnp_bound = 0.0
        for param in params:
            if param < 0.0:
                lnp_bound = -1e20
            # return -np.inf

        # update the size distributions
        dustmodel.set_size_dist(params)

        return dustmodel.lnprob_generic(obsdata) + lnp_bound

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0
        Specific to MCMC fitters (e.g., emcee).

        Parameters
        ----------
        p0 : floats
            Initial values of the parameters
        nwalkers : int
            Number of walkers to initialize

        Returns
        -------
        array of floats
            concatenated set of initial walker positions
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # Initial ball should be in log space
        p = [10**(np.log10(p0) + 1.*np.random.uniform(-1, 1., self.ndim))
             for k in range(self.nwalkers)]

        return p

    def save_results(self, filename, OD, size_dist_uncs=[0]):
        """
        Save fitting results to a file.  Results include the
        size distribution and all predicted observations.

        Creates a FITS file with the results

        Parameters
        ----------
        filename : str
            Name of the file to save the results
        OD : ObsData object
            All the observed data (may not be needed)
        size_dist_uncs : floats
            Uncertainties on the size distributions
        """
        # write a small primary header
        pheader = fits.Header()
        pheader.set('NCOMPS', len(self.components),
                    'number of dust grain components')
        for k, component in enumerate(self.components):
            pheader.set('CNAME'+str(k), component.name,
                        'name of dust grain component')
        pheader.set('SDMODEL', self.sizedisttype,
                    'type of size  distribution')
        pheader.add_comment('Dust Model results written by DustModel.py')
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

            tbhdu = fits.BinTableHDU.from_columns(all_cols)
            tbhdu.header.set('EXTNAME', component.name,
                             'dust grain component name')

            # save the parameter values
            if self.parameters:
                for cparam in self.parameters[component.name].items():
                    tbhdu.header.set(cparam[0], cparam[1],
                                     'parameters of size distribution model')

            hdulist.append(tbhdu)

        # output the resulting observable parameters
        results = self.eff_grain_props(OD, predict_all=True)
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
            results = component.eff_grain_props(OD, predict_all=True)
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

    @staticmethod
    def get_percentile_vals(chain, ndim):
        """
        Compute the 50% +/- 33% values from the samples

        Parameters
        ----------
        chain : sampler.chain
            Chain from the EMCEE sampler
        ndim : int
            number of paramaters

        Returns
        -------
        tuple of floats
            (p50, p84-p50, p50-p16)
        """
        samples = chain.reshape((-1, ndim))
        values = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                     zip(*np.percentile(samples, [16, 50, 84],
                                        axis=0)))
        val_50p, punc, munc = zip(*values)
        return (val_50p, punc, munc)

    def save_50percentile_results(self, oname, sampler, obsdata,
                                  cur_step=None):
        """
        Compute the 50th percentile paramaters, set the size
        distribution, and save the results

        Creates a FITS file with the results

        Parameters
        ----------
        oname : str
            Name of the file to save the results
        sampler : emcee.sampler
            Sampler object from EMCEE run
        obsdata : ObsData object
            All the observed data (may not be needed)
        cur_step : int
            Current step number
        """
        if cur_step is None:
            cur_step = sampler.chain.shape[1]
        fin_size_dist_50p, fin_size_dist_punc, fin_size_dist_munc = \
            self.get_percentile_vals(sampler.chain[:, 0:cur_step+1, :],
                                     self.ndim)
        self.set_size_dist(fin_size_dist_50p)

        # save the model parameters for the size distribution
        # set here so that the saved results have the right info
        self.set_size_dist_parameters(fin_size_dist_50p)

        # save the final size distributions
        self.save_results(oname, obsdata,
                          size_dist_uncs=[fin_size_dist_punc,
                                          fin_size_dist_munc])

    def save_best_results(self, oname, sampler, obsdata,
                          cur_step=None):
        """
        Compute the best fit paramaters, set the size
        distribution, and save the results

        Creates a FITS file with the results

        Parameters
        ----------
        oname : str
            Name of the file to save the results
        sampler : emcee.sampler
            Sampler object from EMCEE run
        obsdata : ObsData object
            All the observed data (may not be needed)
        cur_step : int
            Current step number
        """
        # get the best fit values
        max_lnp = -1e20
        if cur_step is None:
            cur_step = len(sampler.lnprobability[0])
        for k in range(self.nwalkers):
            tmax_lnp = np.max(sampler.lnprobability[k, 0:cur_step])
            if tmax_lnp > max_lnp:
                max_lnp = tmax_lnp
                indxs, = np.where(sampler.lnprobability[k] == tmax_lnp)
                fit_params_best = sampler.chain[k, indxs[0], :]

        self.set_size_dist(fit_params_best)

        # save the best fit size distributions
        self.save_results(oname, obsdata)

# ================================================================


class MRNDustModel(DustModel):
    """
    Dust model that uses powerlaw size distributions with min/max
    sizes (MRN).

    Same keywords and attributes as the parent DustModel class.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sizedisttype = 'MRN'
        self.n_params = [4]*self.n_components
        for component in self.components:
            self.parameters[component.name] = {'C': 1e-25,
                                               'alpha': 3.5,
                                               'a_min': 1e-7,
                                               'a_max': 1e-3}

    def compute_size_dist(self, x, params):
        """
        Compute the size distribution for the input sizes.
        Powerlaw size distribution (aka MRN size distribution)

        sizedist = A*a^-alpha

        where
            a = grain size,
            A = amplitude,
            alpha = exponent of power law,
            amin = min grain size,
            amax = max grain size,

        Parameters
        ----------
        x : floats
            grains sizes
        params : floats
            Size distribution parameters

        Returns
        -------
        floats
            Size distribution as a function of x
        """
        sizedist = params[0]*np.power(x, -1.0*params[1])
        indxs, = np.where(np.logical_or(x < params[2],
                                        x > params[3]))
        if len(indxs) > 0:
            sizedist[indxs] = 0.0

        return sizedist

    def set_size_dist_parameters(self, params):
        """
        Set the size distribution parameters in the object dictonary.

        Parameters
        ----------
        params : floats
            Size distribution parameters
        """
        k1 = 0
        for k, component in enumerate(self.components):
            k2 = k1 + self.n_params[k]
            cparams = params[k1:k2]
            k1 += self.n_params[k]
            self.parameters[component.name] = {'C': cparams[0],
                                               'alpha': cparams[1],
                                               'a_min': cparams[2],
                                               'a_max': cparams[3]}

    @staticmethod
    def lnprob(params, obsdata, dustmodel):
        """
        Compute the ln(prob) given the model parameters

        MRN model paramters for each component are
            A = amplitude
            alpha = negative of the power law exponent
            amin = min grain size
            amax = max grain size

        Parameters
        ----------
        params : array of floats 4 x n_components
            parameters of the MRN model
        obsdata : ObsData object
            observed data for fitting
        dustmodel : DustModel object
            must be passed explicitly as the fitters
            require a static method (is this true?)

        Returns
        -------
        lnprob : float
            natural log of the probability the input parameters
            describe the data
        """
        # priors
        k1 = 0
        lnp_bound = 0.0
        for k, component in enumerate(dustmodel.components):
            # get the parameters for the current component
            k2 = k1 + dustmodel.n_params[k]
            cparams = params[k1:k2]
            k1 += dustmodel.n_params[k]

            # check that amin < amax (params 3 & 4)
            if cparams[2] > cparams[3]:
                lnp_bound = -1e20

            # check that the amin and amax are within the bounds
            # of the dustmodel
            if cparams[2] < component.sizes[0]:
                lnp_bound = -1e20
            if cparams[3] > component.sizes[-1]:
                lnp_bound = -1e20

            # keep the normalization always positive
            if cparams[0] < 0.0:
                lnp_bound = -1e20
            if cparams[1] < 0.0:
                lnp_bound = -1e20

        if lnp_bound < 0.0:
            return lnp_bound
        else:
            dustmodel.set_size_dist(params)
            return dustmodel.lnprob_generic(obsdata) + lnp_bound

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0
        Specific to MCMC fitters (e.g., emcee).

        Parameters
        ----------
        p0 : floats
            Initial values of the parameters
        nwalkers : int
            Number of walkers to initialize

        Returns
        -------
        array of floats
            concatenated set of initial walker positions
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # Initial ball
        # delts = np.array([1.0, 0.01, 1e-8, 1e-5, 1.0, 0.01, 1e-8, 1e-5])
        # p = [p0 + delts*np.random.normal(0., 1., self.ndim)
        p = [10**(np.log10(p0) + 0.1*np.random.uniform(-1, 1., self.ndim))
             for k in range(self.nwalkers)]

        return p

# ================================================================


class WDDustModel(DustModel):
    """
    Dust model that uses the Weingartner & Draine (2001) size distributions.

    Same kewyords and attributes as the parent DustModel class.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sizedisttype = 'WD'

        # set the number of size distribution parametres
        if self.n_components > 0:
            self.n_params = []
            for component in self.components:
                if component.name == 'astro-silicates':
                    self.n_params.append(4)
                    self.parameters['astro-silicates'] = {'C_s': 1.33e-12,
                                                          'a_ts': 0.171e4,
                                                          'alpha_s': -1.41,
                                                          'beta_s': -11.5}
                elif component.name == 'astro-carbonaceous':
                    self.n_params.append(6)
                    self.parameters['astro-carbonaceous'] = {'C_g': 4.15e-11,
                                                             'a_tg': 0.00837e4,
                                                             'alpha_g': -1.91,
                                                             'beta_g': -0.125,
                                                             'a_cg': 0.499e4,
                                                             'b_C': 3.0e-5}
                else:
                    raise ValueError('%s grain material note supported'
                                     % component.name)

    def compute_size_dist(self, x, params):
        """
        Compute the size distribution for the input sizes.

        Parameters
        ----------
        x : floats
            grain sizes
        params : floats
            Size distribution parameters

        Returns
        -------
        floats
            Size distribution as a function of x
        """
        # input grain sizes are in cm, needed in Angstroms
        a = x*1e8

        if len(params) == 6:
            # carbonaceous
            C, a_t, alpha, beta, a_c, input_bC = params
        else:
            # silicates
            C, a_t, alpha, beta = params
            a_c = 0.1e4
            input_bC = None

        # larger grain size distribution
        # same for silicates and carbonaceous grains
        if beta >= 0.0:
            Fa = 1.0 + beta*a/a_t
        else:
            Fa = 1.0/(1.0 - beta*a/a_t)

        Ga = np.full((len(a)), 1.0)
        indxs, = np.where(a > a_t)
        Ga[indxs] = np.exp(-1.0*np.power((a[indxs] - a_t)/a_c, 3.0))

        sizedist = (C/(1e-8*a))*np.power(a/a_t, alpha)*Fa*Ga

        # very small gain size distribution
        # only for carbonaceous grains
        if input_bC is not None:
            a0 = np.array([3.5, 30.])   # in A
            bC = np.array([0.75, 0.25])*input_bC
            sigma = 0.4
            rho = 2.24  # in g/cm^3 for graphite
            mC = 12.0107*1.660e-24

            Da = 0.0
            for i in range(2):
                Bi = ((3.0/(np.power(2.0*np.pi, 1.5)))
                      * (np.exp(-4.5*np.power(sigma, 2.0))
                      / (rho*np.power(1e-8*a0[i], 3.0)*sigma))
                      * (bC[i]*mC
                         / (1.0
                            + erf((3.0*sigma/np.sqrt(2.0))
                                  + np.log(a0[i]/3.5)/(sigma*np.sqrt(2.0))))))

                Da += (Bi/(1e-8*a))*np.exp(-0.5*np.power(np.log(a/a0[i])/sigma,
                                                         2.0))

            sizedist += Da

        return sizedist

    def set_size_dist_parameters(self, params):
        """
        Set the size distribution parameters in the object dictonary.

        Parameters
        ----------
        params : floats
            Size distribution parameters
        """
        k1 = 0
        for k, component in enumerate(self.components):
            k2 = k1 + self.n_params[k]
            cparams = params[k1:k2]
            k1 += self.n_params[k]
            if component.name == 'astro-silicates':
                self.parameters['astro-silicates'] = {'C_s': cparams[0],
                                                      'a_ts': cparams[1],
                                                      'alpha_s': cparams[2],
                                                      'beta_s': cparams[3]}
            else:
                self.parameters['astro-carbonaceous'] = {'C_g': cparams[0],
                                                         'a_tg': cparams[1],
                                                         'alpha_g': cparams[2],
                                                         'beta_g': cparams[3],
                                                         'a_cg': cparams[4],
                                                         'b_C': cparams[5]}

    @staticmethod
    def lnprob(params, obsdata, dustmodel):
        """
        Compute the ln(prob) given the model parameters

        Parameters
        ----------
        params : array of floats 4
            parameters of the WD model
        obsdata : ObsData object
            observed data for fitting
        dustmodel : DustModel object
            must be passed explicitly as the fitters
            require a static method (is this true?)

        Returns
        -------
        lnprob : float
            natural log of the probability the input parameters
            describe the data
        """
        # priors
        k1 = 0
        lnp_bound = 0.0
        for k, component in enumerate(dustmodel.components):
            # get the parameters for the current component
            k2 = k1 + dustmodel.n_params[k]
            cparams = params[k1:k2]
            k1 += dustmodel.n_params[k]

            # keep the normalization always positive
            if cparams[0] < 0.0:
                lnp_bound = -1e20
            if cparams[1] < 0.0:
                lnp_bound = -1e20
            if len(cparams) == 6:
                if cparams[4] < 0.0:
                    lnp_bound = -1e20
                if cparams[5] < 0.0:
                    lnp_bound = -1e20

        if lnp_bound < 0.0:
            return lnp_bound
        else:
            dustmodel.set_size_dist(params)
            return dustmodel.lnprob_generic(obsdata) + lnp_bound

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0
        Specific to MCMC fitters (e.g., emcee).

        Parameters
        ----------
        p0 : floats
            Initial values of the parameters
        nwalkers : int
            Number of walkers to initialize

        Returns
        -------
        array of floats
            concatenated set of initial walker positions
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # some parameters are negative, so need to be handled
        psigns = np.sign(p0)
        p = [psigns*(10**(np.log10(np.absolute(p0))
                     + 0.1*np.random.uniform(-1, 1., self.ndim)))
             for k in range(self.nwalkers)]

        return p
