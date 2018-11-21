from __future__ import print_function

import math
import numpy as np


__all__ = ['lnprob_all', 'get_percentile_vals',
           'mrn_size_model']


def lnprob_all(obsdata, dustmodel):
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


def get_percentile_vals(chain, ndim):
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


def mrn_size_model(a, params):
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
