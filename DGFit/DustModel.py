import numpy as np
from astropy.io import fits

from DGFit.DustGrains import DustGrains
from DGFit.DMhelpers import (lnprob_all,
                             get_percentile_vals,
                             mrn_size_model)

__all__ = ['DustModel',
           'MRNDustModel']


class DustModel():
    """
    Full dust model including arbirary size and composition distributions.
    Includes the physical properties of the individual dust grains.

    Dust model that has each bin as an independent variable in the
    grain size distribution providing a truly arbitrary specification.

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
    """
    def __init__(self):
        self.origin = None
        self.n_components = 0
        self.components = []
        self.sizedisttype = 'bins'

    @staticmethod
    def lnprob(params, obsdata, dustmodel):

        # make sure the size distributions are all positve
        lnp_bound = 0.0
        for param in params:
            if param < 0.0:
                lnp_bound = -1e20
            # return -np.inf

        # update the size distributions
        #  the input params are the concatenated size distributions
        dustmodel.set_size_dist(params)

        return lnprob_all(obsdata, dustmodel) + lnp_bound

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # Initial ball should be in log space
        p = [10**(np.log10(p0) + 1.*np.random.uniform(-1, 1., self.ndim))
             for k in range(self.nwalkers)]

        # ensure that all the walkers start with positive values
        for pc in p:
            for pcs in pc:
                if pcs <= 0.0:
                    pcs = 1e-20

        return p

    def save_percentile_vals(self, oname, sampler, obsdata,
                             cur_step=None):
        """
        Save the percentile values using the sampler chain
        """
        if cur_step is None:
            cur_step = sampler.chain.shape[1]
        fin_size_dist_50p, fin_size_dist_punc, fin_size_dist_munc = \
            get_percentile_vals(sampler.chain[:, 0:cur_step+1, :], self.ndim)
        self.set_size_dist(fin_size_dist_50p)

        # save the final size distributions
        self.save(oname, obsdata,
                  size_dist_uncs=[fin_size_dist_punc, fin_size_dist_munc])

    def save_best_vals(self, oname, sampler, obsdata,
                       cur_step=None):
        """
        Save the best fit values using the sampler chain
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
        self.save(oname, obsdata)

    def predict_full_grid(self, componentnames, path='./'):
        """
        Read in the precomputed dust grain physical properties from files.

        Parameters
        ----------
        componentnames : list of strings
            names of dust grain materials
        path : type
            path to files

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
                              path=path)
            self.components.append(cur_DG)

    def predict_on_observed_data(self, full_dustmodel, observeddata):
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


class MRNDustModel(DustModel):
    """
    Dust model that uses powerlaw size distributions with min/max
    sizes (MRN).
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)
        self.sizedisttype = 'MRN'

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
        k1 = 0
        n_mrn_params = 4
        lnp_bound = 0.0
        for component in dustmodel.components:
            k2 = k1 + n_mrn_params
            cparams = params[k1:k2]
            component.size_dist[:] = mrn_size_model(component.sizes[:],
                                                    cparams)
            k1 += n_mrn_params
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
            return lnprob_all(obsdata, dustmodel) + lnp_bound

    def set_size_dist(self, params):
        """
        Set the size distributions of the components based
        on the input parameters

        Parameters
        ----------
        params : array of floats
            MRN parameters (4 per component)

        Returns
        -------
        updated size_distributions of the componenets
        """
        k1 = 0
        n_mrn_params = 4
        for component in self.components:
            k2 = k1 + n_mrn_params
            cparams = params[k1:k2]
            component.size_dist[:] = mrn_size_model(component.sizes[:],
                                                    cparams)
            k1 += n_mrn_params

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0.
        Specific to MCMC fitters (e.g., emcee).

        Parameters
        ----------
        p0 : type
            Description of parameter `p0`.
        nwalkers : type
            Description of parameter `nwalkers`.

        Returns
        -------
        type
            Description of returned object.
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # Initial ball
        # delts = np.array([1.0, 0.01, 1e-8, 1e-5, 1.0, 0.01, 1e-8, 1e-5])
        # p = [p0 + delts*np.random.normal(0., 1., self.ndim)
        p = [10**(np.log10(p0) + 0.1*np.random.uniform(-1, 1., self.ndim))
             for k in range(self.nwalkers)]

        return p

    def save_percentile_vals(self, oname, sampler, obsdata,
                             cur_step=None):
        """
        Save the percentile values using the sampler chain

        Parameters
        ----------
        oname : type
            Description of parameter `oname`.
        sampler : type
            Description of parameter `sampler`.
        obsdata : type
            Description of parameter `obsdata`.
        cur_step : type
            Description of parameter `cur_step`.

        Returns
        -------
        type
            Description of returned object.

        """
        if cur_step is None:
            cur_step = sampler.chain.shape[1]
        fin_param_50p, fin_param_punc, fin_param_munc = \
            get_percentile_vals(sampler.chain[:, 0:cur_step+1, :], self.ndim)
        self.set_size_dist(fin_param_50p)

        # save the final size distributions
        self.save(oname, obsdata)

    def save_best_vals(self, oname, sampler, obsdata,
                       cur_step=None):
        """
        Save the best fit values using the sampler chain

        Parameters
        ----------
        oname : type
            Description of parameter `oname`.
        sampler : type
            Description of parameter `sampler`.
        obsdata : type
            Description of parameter `obsdata`.
        cur_step : type
            Description of parameter `cur_step`.

        Returns
        -------
        type
            Description of returned object.
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
        self.save(oname, obsdata)
