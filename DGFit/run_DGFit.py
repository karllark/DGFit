#!/usr/bin/env python
#
# DGFit - program to determine dust grain size & composition from
#         fitting dust grain observables
#  dust observables = extinction, scattering properties, depletions,
#                     emission, etc.
#
# Started: Jan 2015 (KDG)
#  added IR emission: Feb 2015 (KDG)
#  udpated to move plotting to a separate function: Dec 2015 (KDG)
#
from __future__ import print_function

import sys
import math
import time
import argparse

import numpy as np
from astropy.io import fits

# from scipy.interpolate import interp1d
from scipy.optimize import minimize

# from lmfit import minimize, Parameters

import emcee

from DGFit.DustModel import DustModel
from DGFit.ObsData import ObsData
from DGFit.drainesizedist import WGsizedist_graphite
# import ObsData_Azv18 as ObsData


# get the ln(prob) for the dust grain size/composition distribution
#  defined in dustmodel
def lnprob_all(obsdata, dustmodel):

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
            # hard limit at 1.5x the total possible abundaces
            #      (all atoms in dust)
            # if natoms[atomname] > 1.5*obsdata.total_abundance[atomname][0]:
                # print('boundary issue')
                # return -np.inf
                # pass
            # only add if natoms > depletions
            # elif natoms[atomname] > obsdata.abundance[atomname][0]:
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


# compute the ln(prob) for discrete size distributions
def lnprob_discrete(params, obsdata, dustmodel):

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


# compute the ln(prob) for discrete size distributions
def lnprob_WGsizedist(params, obsdata, dustmodel):

    # update the size distributions
    #   params are the WG size dist parameters
    # params = bC_input, Cg, a_tg, a_cg, alpha_g, beta_g
    for component in dustmodel.components:
        csizedist = WGsizedist_graphite(component.sizes, params)
    dustmodel.set_size_dist_WG(csizedist)

    return lnprob_all(obsdata, dustmodel)


def DGFit_cmdparser():

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fast",
                        help="Use minimal walkers, steps, burns to debug code",
                        action="store_true")
    parser.add_argument("-s", "--slow",
                        help="Use lots of walkers, n_steps, n_burn",
                        action="store_true")
    parser.add_argument("--limit_abund", action="store_true",
                        help="Limit based on abundances")
    parser.add_argument("--usemin", action="store_true",
                        help="Find min before EMCEE")
    parser.add_argument("-r", "--read", default="",
                        help="Read size distribution from disk")
    parser.add_argument("-t", "--tag", default='dgfit_test',
                        help="basename to use for output files")
    parser.add_argument("-c", "--cpus", metavar=int, default=4,
                        help="number of cpus to use")
    parser.add_argument("--nolarge", action="store_true",
                        help="Deweight a > 0.5 micron by 1e-10")
    parser.add_argument("--smc", help="use an SMC sightline",
                        action="store_true")

    return parser


# main fitting code
if __name__ == "__main__":

    parser = DGFit_cmdparser()

    args = parser.parse_args()

    # set the basename of the output
    basename = args.tag

    # save the start time
    start_time = time.clock()

    # get the observed data
    if args.smc:
        path = 'DGFit/data/smc_azv215'
        obsdata = ObsData('%s/azv215_50p_ext.fits' % path,
                          '%s/azv215_avnhi.dat' % path,
                          '%s/SMC_AzV215_abundances.dat' % path,
                          None,
                          None)
    else:
        path = 'DGFit/data/mw_rv31'
        obsdata = ObsData(['%s/MW_diffuse_Gordon09_band_ext.dat' % path,
                           '%s/MW_diffuse_Gordon09_iue_ext.dat' % path,
                           '%s/MW_diffuse_Gordon09_fuse_ext.dat' % path],
                          '%s/MW_diffuse_Gordon09_avnhi.dat' % path,
                          '%s/MW_diffuse_Jenkins09_abundances.dat' % path,
                          '%s/MW_diffuse_Compiegne11_ir_emission.dat' % path,
                          '%s/dust_scat.dat' % path,
                          ext_tags=['band', 'iue', 'fuse'],
                          scat_path='%s/Scat_Data/' % path)

    # get the dust model on the full wavelength grid
    dustmodel_full = DustModel()
    dustmodel_full.predict_full_grid(['astro-silicates', 'astro-carbonaceous'],
                                     path='DGFit/data/indiv_grain/')

    # get the dust model predicted on the observed data grids
    dustmodel = DustModel()
    dustmodel.predict_observed_data(dustmodel_full, obsdata)

    # replace the default size distribution with one from a file
    if args.read != "":
        for k, component in enumerate(dustmodel.components):
            fitsdata = fits.getdata(args.read, k+1)
            if len(component.size_dist) != len(fitsdata[:][1]):
                component.size_dist = 10**np.interp(np.log10(component.sizes),
                                                    np.log10(fitsdata['SIZE']),
                                                    np.log10(fitsdata['DIST']))
            else:
                component.size_dist = fitsdata['DIST']

    else:
        # check that the default size distributions give approximately
        #     the right level of the A(lambda)/N(HI) curve
        # if not, adjust the overall level of the size distributions to
        #     get them close
        results = dustmodel.eff_grain_props(obsdata)
        cabs = results['cabs']
        csca = results['csca']
        dust_alnhi = 1.086*(cabs + csca)
        ave_model = np.average(dust_alnhi)
        ave_data = np.average(obsdata.ext_alnhi)
        ave_ratio = ave_data/ave_model
        if (ave_ratio < 0.5) | (ave_ratio > 2):
            for component in dustmodel.components:
                component.size_dist *= ave_ratio

        # results = dustmodel.eff_grain_props()
        # natoms = results[2]
        # max_violation = 0.0
        # for atomname in natoms.keys():
            # cur_violation = natoms[atomname]/obsdata.abundance[atomname][0]
            # if cur_violation > max_violation:
            #    max_violation = cur_violation

        # if max_violation > 2:
        #    for component in dustmodel.components:
        #        component.size_dist *= 1.9/max_violation

            # deweight large grains (test)
    if args.nolarge:
        indxs, = np.where(component.sizes > 0.5e-4)
        if len(indxs) > 0:
            print('deweighting sizes > 0.5 micron')
            component.size_dist[indxs] *= 1e-10

    # save the starting model
    dustmodel.save(basename + '_sizedist_start.fits', obsdata)

    # setup time
    setup_time = time.clock()
    print('setup time taken: ', (setup_time - start_time)/60., ' min')

    # inital guesses at parameters
    p0 = dustmodel.components[0].size_dist
    for k in range(1, dustmodel.n_components):
        p0 = np.concatenate([p0, dustmodel.components[k].size_dist])

    # call scipy.optimize to get a better initial guess
    if args.usemin:
        print(p0)
        print(lnprob_discrete(p0, obsdata, dustmodel))

        # generate the bounds
        p0_bounds = []
        for k in range(len(p0)):
            p0_bounds.append((0.0, 1e20))

            # setup what can be fit
        obsdata.fit_extinction = True
        obsdata.fit_abundance = False
        obsdata.fit_ir_emission = False
        obsdata.fit_scat_a = False
        obsdata.fit_scat_g = False

        # neg_lnprobsed = lambda *args: -1.0*lnprob_discrete(*args)
        def neg_lnprobsed(*args): -1.0*lnprob_discrete(*args)
        better_start = minimize(neg_lnprobsed, p0, args=(obsdata, dustmodel),
                                bounds=p0_bounds, method='L-BFGS-B')
        print(better_start.success)
        print(better_start.x)
        exit()

    # import scipy.optimize as op
    # nll = lambda *args: -lnprobsed(*args)
    # result = op.minimize(nll, p0, args=(obsdata, dustmodel))
    # print(result)
    # exit()

    # trying with lmfit (not tested)
    # params = Parameters()
    # for i in range(n_params):
    #   params.add('p'+str(i),value=p0[i],min=0.0)
    # out = minimize(kext_residuals, params, args=(1.0/xdata, ydata))
    # print(out.params)

    ndim = len(p0)
    print('# params = ', ndim)
    if args.fast:
        print('using the fast params')
        nwalkers = 2*ndim
        nsteps = 100
        burn = 50
    elif args.slow:
        print('using the slow params')
        nwalkers = 2*ndim
        nsteps = 10000
        burn = 5000
    else:
        nwalkers = 2*ndim
        nsteps = 1000
        burn = 500

    # setting up the walkers to start "near" the inital guess
    #   and initial ball should be in log space
    p = [10**(np.log10(p0) + 1.*np.random.uniform(-1, 1., ndim))
         for k in range(nwalkers)]

    # ensure that all the walkers start with positive values
    for pc in p:
        for pcs in pc:
            if pcs <= 0.0:
                pcs = 0.0

    # make sure each walker starts with allowed abundances
    if args.limit_abund:
        for pc in p:
            dustmodel.set_size_dist(pc)
            results = dustmodel.eff_grain_props(obsdata)
            cabs = results['cabs']
            csca = results['csca']
            natoms = results['natoms']
            max_violation = 0.0
            for atomname in natoms.keys():
                cur_violation = (natoms[atomname]
                                 / (obsdata.abundance[atomname][0] +
                                    obsdata.abundance[atomname][1]))
                if cur_violation > max_violation:
                    max_violation = cur_violation
            if max_violation > 2:
                pc *= 1.9/max_violation

    # setup the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_discrete,
                                    args=(obsdata, dustmodel),
                                    threads=int(args.cpus))

    # burn in the walkers
    # pos, prob, state = sampler.run_mcmc(p, burn)
    print("burn")
    width = 60
    for i, result in enumerate(sampler.sample(p, iterations=burn)):
        n = int((width+1) * float(i) / burn)
        sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
    sys.stdout.write("\n")

    pos, prob, state = result

    # rest the sampler
    sampler.reset()

    # do the full sampling
    # pos, prob, state = sampler.run_mcmc(pos, nsteps, rstate0=state)

    print("afterburn")
    width = 60
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps,
                                              rstate0=state)):
        n = int((width+1) * float(i) / nsteps)
        sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
    sys.stdout.write("\n")

    pos, prob, state = result

    # untested code from emcee webpages for incrementally saving the chains
    # f = open("chain.dat", "w")
    # f.close()
    # for k, result in enumerate(sampler.sample(pos0, iterations=500,
    #                           storechain=False)):
    #    print(k)
    #    position = result[0]
    #    f = open("chain.dat", "a")
    #    for k in range(position.shape[0]):
    #        f.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
    #    f.close()

    emcee_time = time.clock()
    print('emcee time taken: ', (emcee_time - setup_time)/60., ' min')

    # get the best fit values
    max_lnp = -1e20
    for k in range(nwalkers):
        tmax_lnp = np.max(sampler.lnprobability[k])
        if tmax_lnp > max_lnp:
            max_lnp = tmax_lnp
            indxs, = np.where(sampler.lnprobability[k] == tmax_lnp)
            fit_params_best = sampler.chain[k, indxs[0], :]

    dustmodel.set_size_dist(fit_params_best)

    # save the best fit size distributions
    dustmodel.save(basename + '_sizedist_best.fits', obsdata)

    # get the 50p values and uncertainties
    samples = sampler.chain.reshape((-1, ndim))
    values = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                 zip(*np.percentile(samples, [16, 50, 84],
                                    axis=0)))
    fin_size_dist_50p, fin_size_dist_punc, fin_size_dist_munc = zip(*values)

    # 50p dust params
    dustmodel.set_size_dist(fin_size_dist_50p)

    # save the final size distributions
    dustmodel.save(basename + '_sizedist.fits', obsdata,
                   size_dist_uncs=[fin_size_dist_punc, fin_size_dist_munc])
