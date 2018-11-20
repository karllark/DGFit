from __future__ import print_function

import math
import numpy as np


__all__ = ['DGFit_bins', 'DGFit_MRN']


def _lnprob_all(obsdata, dustmodel):
    """
    Compute the ln(prob) for the dust grain size and composition
    distribution as defined in the dustmodel
    """
    # get the integrated dust properties
    results = dustmodel.eff_grain_props(obsdata)

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


def _get_percentile_vals(chain, ndim):
    """
    Compute the 50% +/- 33% values from the samples
    """
    # get the 50p values and uncertainties
    samples = chain.reshape((-1, ndim))
    values = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                 zip(*np.percentile(samples, [16, 50, 84],
                                    axis=0)))
    val_50p, punc, munc = zip(*values)
    return (val_50p, punc, munc)


class DGFit_bins():
    """
    Model that parametrizes the dust size distributions as discrete bins
    """
    def __init__(self):
        self.type = 'bins'

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

        return _lnprob_all(obsdata, dustmodel) + lnp_bound

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

    def save_percentile_vals(self, oname, dustmodel, sampler, obsdata,
                             cur_step=None):
        """
        Save the percentile values using the sampler chain
        """
        if cur_step is None:
            cur_step = sampler.chain.shape[1]
        fin_size_dist_50p, fin_size_dist_punc, fin_size_dist_munc = \
            _get_percentile_vals(sampler.chain[:, 0:cur_step+1, :], self.ndim)
        dustmodel.set_size_dist(fin_size_dist_50p)

        # save the final size distributions
        dustmodel.save(oname, obsdata,
                       size_dist_uncs=[fin_size_dist_punc, fin_size_dist_munc])

    def save_best_vals(self, oname, dustmodel, sampler, obsdata,
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

        dustmodel.set_size_dist(fit_params_best)

        # save the best fit size distributions
        dustmodel.save(oname, obsdata)


def _mrn_size_model(a, params):
    """
    MRN size distribution

    MRN model paramters are
        A = amplitude
        alpha = exponent of power law
        amin = min grain size
        amax = max grain size

    sizedist = A*a^-alpha
    """
    sizedist = params[0]*np.power(a, -1.0*params[1])
    indxs, = np.where(np.logical_or(a < params[2],
                                    a > params[3]))
    if len(indxs) > 0:
        sizedist[indxs] = 0.0

    return(sizedist)


class DGFit_MRN():
    """
    Model that uses the Weingartner & Draine functional form for the
    size distributions
    """
    def __init__(self):
        self.type = 'MRN'

    @staticmethod
    def lnprob(params, obsdata, dustmodel):
        """
        Compute the ln(prob) given the model parameters

        MRN model paramters for each component are
            A = amplitude
            alpha = negative of the power law exponent
            amin = min grain size
            amax = max grain size
        """
        k1 = 0
        n_mrn_params = 4
        lnp_bound = 0.0
        for component in dustmodel.components:
            k2 = k1 + n_mrn_params
            cparams = params[k1:k2]
            component.size_dist[:] = _mrn_size_model(component.sizes[:],
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
            return _lnprob_all(obsdata, dustmodel) + lnp_bound

    def set_size_dist(self, params, dustmodel):
        k1 = 0
        n_mrn_params = 4
        for component in dustmodel.components:
            k2 = k1 + n_mrn_params
            cparams = params[k1:k2]
            component.size_dist[:] = _mrn_size_model(component.sizes[:],
                                                     cparams)
            k1 += n_mrn_params

    def initial_walkers(self, p0, nwalkers):
        """
        Setup the walkers based on the initial parameters p0
        """
        self.ndim = len(p0)
        self.nwalkers = nwalkers
        # Initial ball
        # delts = np.array([1.0, 0.01, 1e-8, 1e-5, 1.0, 0.01, 1e-8, 1e-5])
        # p = [p0 + delts*np.random.normal(0., 1., self.ndim)
        p = [10**(np.log10(p0) + 0.1*np.random.uniform(-1, 1., self.ndim))
             for k in range(self.nwalkers)]

        return p

    def save_percentile_vals(self, oname, dustmodel, sampler, obsdata,
                             cur_step=None):
        """
        Save the percentile values using the sampler chain
        """
        if cur_step is None:
            cur_step = sampler.chain.shape[1]
        fin_param_50p, fin_param_punc, fin_param_munc = \
            _get_percentile_vals(sampler.chain[:, 0:cur_step+1, :], self.ndim)
        self.set_size_dist(fin_param_50p, dustmodel)

        # save the final size distributions
        dustmodel.save(oname, obsdata)

    def save_best_vals(self, oname, dustmodel, sampler, obsdata,
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

        self.set_size_dist(fit_params_best, dustmodel)

        # save the best fit size distributions
        dustmodel.save(oname, obsdata)


class DGFit_WG():
    """
    Model that uses the Weingartner & Draine functional form for the
    size distributions
    """
    def __init__(self):
        self.type = 'WG'

    def lnprob(params, obsdata, dustmodel):
        # update the size distributions
        #   params are the WG size dist parameters
        # params = bC_input, Cg, a_tg, a_cg, alpha_g, beta_g
        for component in dustmodel.components:
            pass
            # csizedist = WGsizedist_graphite(component.sizes, params)
            # dustmodel.set_size_dist_WG(csizedist)

        return _lnprob_all(obsdata, dustmodel)
