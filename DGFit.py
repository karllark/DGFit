#!/usr/bin/env python2.7
#
# DGFit - program to determine dust grain size & composition from fitting dust grain observables
#  dust observables = extinction, scattering properties, depletions, emission, etc.
#
# Started: Jan 2015 (KDG)
#  added IR emission: Feb 2015 (KDG)
# 
from __future__ import print_function

import math
import time
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
from astropy.io import fits

import emcee
import triangle

import DustModel
import ObsData

# compute the ln(prob) for an input set of model parameters
def lnprobsed(params, ObsData, DustModel):

    # make sure the size distributions are all positve
    for param in params:
        if param < 0.0:
            return -np.inf

    # update the size distributions
    #  the input params are the concatenated size distributions
    DustModel.set_size_dist(params)

    #print(params)

    # get the integrated dust properties
    cabs, csca, natoms, emission = DustModel.eff_grain_props()
    cext = cabs + csca
    dust_alnhi = 1.086*cext
    
    # compute the ln(prob) for A(l)/N(HI)
    lnp_alnhi = -0.5*np.sum((((obsdata.alnhi - dust_alnhi)/(0.10*obsdata.alnhi))**2))
    #lnp_alnhi /= obsdata.n_wavelengths

    # compute the ln(prob) for the depletions
    lnp_dep = 0.0
    for atomname in natoms.keys():
        if natoms[atomname] > 2.*obsdata.depletions[atomname][0]: # hard limit at 2x depletions
            return -np.inf
        elif natoms[atomname] > obsdata.depletions[atomname][0]: # only add if natoms > depletions
            lnp_dep = ((natoms[atomname] - obsdata.depletions[atomname][0])/obsdata.depletions[atomname][1])**2
    lnp_dep *= -0.5
    #lnp_dep /= len(natoms.keys())

    # compute the ln(prob) for IR emission
    lnp_emission = -0.5*np.sum((((obsdata.ir_emission - emission[obsdata.ir_emission_indxs])/(obsdata.ir_emission_uncs))**2))

    # compute the ln(prob) for the dust albedo
    albedo = csca/cext
    lnp_albedo = -0.5*np.sum((((obsdata.scat_albedo - albedo[obsdata.scat_indxs])/(obsdata.scat_albedo_unc))**2))

    # combine the lnps
    lnp = lnp_alnhi + lnp_dep + lnp_emission + lnp_albedo
    
    if math.isinf(lnp) | math.isnan(lnp):
        print(lnp)
        print(params)
        exit()
    else:
        return lnp

# main fitting code
if __name__ == "__main__":
    
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fast", help="Use minimal walkers, n_steps, n_burn to debug code",
                        action="store_true")
    parser.add_argument("-s", "--slow", help="Use lots of walkers, n_steps, n_burn",
                        action="store_true")
    parser.add_argument("-r", "--read", help="Read size distribution from disk",
                        action="store_true")
    parser.add_argument("-t", "--triangle", help="Plot the 'triangle' plot showing the 1D and 2D likelihood functions",
                        action="store_true")
    parser.add_argument("-w", "--walkers", help="Plot the walker values",
                        action="store_true")
    parser.add_argument("-o", "--plot_old", help="Plot the starting values from the original starting ones",
                        action="store_true")
    args = parser.parse_args()

    # save the start time 
    start_time = time.clock()

    # get the dust model 
    min_wave = 0.09
    max_wave = 3.0
    min_wave_emission = 1.0
    max_wave_emission = 1000.0
    dustmodel = DustModel.DustModel(['astro-silicates','astro-carbonaceous'], path='/home/kgordon/Dirty_v2/write_grain/indiv_grain/',
                                    min_wave=min_wave,max_wave=max_wave,
                                    min_wave_emission=min_wave_emission,max_wave_emission=max_wave_emission)

    # get the observed data to fit
    obsdata = ObsData.ObsData(dustmodel.components[0].wavelengths, dustmodel.components[0].wavelengths_emission)

    if args.plot_old:    
        # check that the default size distributions give approximately the right level of the A(lambda)/N(HI) curve
        #  if not, adjust the overall level of the size distributions to get them close
        cabs, csca, natoms, emission = dustmodel.eff_grain_props()
        dust_alnhi = 1.086*(cabs + csca)
        ave_model = np.average(dust_alnhi)
        ave_data = np.average(obsdata.alnhi)
        ave_ratio = ave_data/ave_model
        if (ave_ratio < 0.5) | (ave_ratio > 2):
            for component in dustmodel.components:
                component.size_dist *= ave_ratio

        # temp fix for carbon starting way too high
        dustmodel.components[0].size_dist *= 1.75
        dustmodel.components[1].size_dist /= 4.

        # inital guesses at parameters
        p0_start = dustmodel.components[0].size_dist
        for k in range(1,dustmodel.n_components):
            p0_start = np.concatenate([p0_start,dustmodel.components[k].size_dist])

        # save the starting model curve
        cabs_start, csca_start, natoms_start, emission_start = dustmodel.eff_grain_props()
        dust_alnhi_start = 1.086*(cabs_start + csca_start)

    if args.read:
        for k, component in enumerate(dustmodel.components):
            fitsdata = fits.getdata('dgfit_test_sizedist.fits',k+1)
            for i in range(component.n_sizes):
                component.size_dist[i] = fitsdata[i][1]
                #if component.sizes[i] > 0.5e-4:
                #    component.size_dist[i] *= 1e-4
    else:
        # check that the default size distributions give approximately the right level of the A(lambda)/N(HI) curve
        #  if not, adjust the overall level of the size distributions to get them close
        cabs, csca, natoms, emission = dustmodel.eff_grain_props()
        dust_alnhi = 1.086*(cabs + csca)
        ave_model = np.average(dust_alnhi)
        ave_data = np.average(obsdata.alnhi)
        ave_ratio = ave_data/ave_model
        if (ave_ratio < 0.5) | (ave_ratio > 2):
            for component in dustmodel.components:
                component.size_dist *= ave_ratio

        cabs, csca, natoms, emission = dustmodel.eff_grain_props()
        max_violation = 0.0
        for atomname in natoms.keys():
            cur_violation = natoms[atomname]/obsdata.depletions[atomname][0]
            if cur_violation > max_violation:
                max_violation = cur_violation
        if max_violation > 2:
            for component in dustmodel.components:
                component.size_dist *= 1.9/max_violation        

    if not args.plot_old:
        # save the starting model curve
        cabs_start, csca_start, natoms_start, emission_start = dustmodel.eff_grain_props()
        dust_alnhi_start = 1.086*(cabs_start + csca_start)
    
    # setup time
    setup_time = time.clock()
    print('setup time taken: ',(setup_time - start_time)/60., ' min')

    # inital guesses at parameters
    p0 = dustmodel.components[0].size_dist
    for k in range(1,dustmodel.n_components):
        p0 = np.concatenate([p0,dustmodel.components[k].size_dist])

    if not args.plot_old:
        p0_start = p0

    ndim = len(p0)
    print('# params = ', ndim)
    if args.fast:
        print('using the fast params')
        nwalkers = 2*ndim
        nsteps = 50
        burn   = 5
    elif args.slow:
        print('using the slow params')
        nwalkers = 2*ndim
        nsteps = 5000
        burn   = 200000
    else:
        nwalkers = 2*ndim
        nsteps = 1000
        burn   = 1000

    # setting up the walkers to start "near" the inital guess
    p  = [ p0*(1+0.25*np.random.normal(0,1.,ndim))  for k in xrange(nwalkers)]

    # ensure that all the walkers start with positive values
    for pc in p:
        for pcs in pc:
            if pcs <= 0.0:
                pcs = 0.0

    # make sure each walker starts with allowed abundances
    for pc in p:
        dustmodel.set_size_dist(pc)
        cabs, csca, natoms, emission = dustmodel.eff_grain_props()
        max_violation = 0.0
        for atomname in natoms.keys():
            cur_violation = natoms[atomname]/obsdata.depletions[atomname][0]
            if cur_violation > max_violation:
                max_violation = cur_violation
        if max_violation > 2:
            pc *= 1.9/max_violation        
            
    # setup the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobsed, args=(obsdata, dustmodel), threads=4)

    # burn in the walkers
    pos, prob, state = sampler.run_mcmc(p, burn)

    # rest the sampler
    sampler.reset()

    # do the full sampling
    pos, prob, state = sampler.run_mcmc(pos, nsteps, rstate0=state)

    emcee_time = time.clock()
    print('emcee time taken: ',(emcee_time - setup_time)/60., ' min')

    # get the best fit values
    max_lnp = -1e6
    for k in range(nwalkers):
        tmax_lnp = np.max(sampler.lnprobability[k])
        if tmax_lnp > max_lnp:
            max_lnp = tmax_lnp
            indxs, = np.where(sampler.lnprobability[k] == tmax_lnp)
            fit_params_best = sampler.chain[k,indxs[0],:]

    dustmodel.set_size_dist(fit_params_best)
    fin_size_dist_best = fit_params_best
    cabs_best, csca_best, natoms_best, emission_best = dustmodel.eff_grain_props()
    dust_alnhi_best = 1.086*(cabs_best + csca_best)

    # get the 50p values and uncertainties
    samples = sampler.chain.reshape((-1, ndim))
    values = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                 zip(*np.percentile(samples, [16, 50, 84],
                                    axis=0)))
    fin_size_dist_50p, fin_size_dist_punc, fin_size_dist_munc = zip(*values)

    # 50p dust params
    dustmodel.set_size_dist(fin_size_dist_50p)
    cabs, csca, natoms, emission = dustmodel.eff_grain_props()
    dust_alnhi = 1.086*(cabs + csca)

    # plot the observed, initial and final extinction curves
    fig, ax = pyplot.subplots(nrows=2,ncols=3, sharex=False, figsize=(20,10))

    ax[0,0].plot(obsdata.wavelengths, obsdata.alnhi, 'k-', label='Observed')
    ax[0,0].plot(obsdata.wavelengths, dust_alnhi, 'b-', label='50p fit parameters')
    ax[0,0].plot(obsdata.wavelengths, dust_alnhi_best, 'g-', label='best fit parameters')
    ax[0,0].plot(obsdata.wavelengths, dust_alnhi_start, 'r--', label='starting parameters')
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    ax[0,0].set_xlabel('$\lambda$ [$\mu$m]')
    ax[0,0].set_ylabel('$A(\lambda)/N(HI)$')
    ax[0,0].set_xlim(0.9*min_wave,1.1*max_wave)
    ax[0,0].legend()

    # plot the albedos
    ax[0,2].plot(obsdata.wavelengths, csca/(cabs+csca), 'b-', label='50p fit parameters')
    ax[0,2].plot(obsdata.wavelengths, csca_best/(cabs_best+csca_best), 'g-', label='best fit parameters')
    ax[0,2].plot(obsdata.wavelengths, csca_start/(cabs_start+csca_start), 'r--', label='starting parameters')
    ax[0,2].errorbar(obsdata.scat_waves, obsdata.scat_albedo, yerr=obsdata.scat_albedo_unc, fmt='ko', label='Observed')
    ax[0,2].set_xscale('log')
    ax[0,2].set_xlabel('$\lambda$ [$\mu$m]')
    ax[0,2].set_ylabel('albedo')
    ax[0,2].set_xlim(0.9*min_wave,1.1*max_wave)

    # plot the emission
    emis_waves = dustmodel.components[0].wavelengths_emission
    ax[0,1].plot(emis_waves, emission, 'b-', label='50p fit parameters')
    ax[0,1].plot(emis_waves, emission_best, 'g-', label='best fit parameters')
    ax[0,1].plot(emis_waves, emission_start, 'r--', label='starting parameters')
    ax[0,1].plot(obsdata.ir_emission_waves, obsdata.ir_emission, 'ko', label='Observed')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlabel('$\lambda$ [$\mu$m]')
    ax[0,1].set_ylabel('Emission')
    ax[0,1].set_xlim(0.9*min_wave_emission,1.1*max_wave_emission)

    # plot the atomic abundances used
    n_atoms = len(natoms.values())
    aindxs = np.arange(n_atoms)
    width = 0.35
    ax[1,0].bar(aindxs, natoms_start.values(), width, color='r')
    ax[1,0].bar(aindxs+0.5*width, natoms_best.values(), width, color='g')
    ax[1,0].bar(aindxs+1.0*width, natoms.values(), width, color='b')

    ax[1,0].errorbar(aindxs+width, [obsdata.depletions[x][0] for x in natoms.keys()],
                     yerr=[obsdata.depletions[x][1] for x in natoms.keys()], fmt='ko')

    ax[1,0].set_ylabel('$N(X)/N(HI)$')
    ax[1,0].set_xticks(aindxs+(0.75*width))
    ax[1,0].set_xticklabels( natoms.keys() )
    
    # plot the final dust size distributions
    k1 = 0
    pfmts = ['b--','r--','g--']
    pfmts_fin = ['b-','r-','g-']
    pfmts_best = ['bo','ro','go']
    for i, component in enumerate(dustmodel.components):
        k2 = k1 + component.n_sizes
        sizes = component.sizes
        size_dist_start = p0_start[k1:k2]
        size_dist = fin_size_dist_50p[k1:k2]
        size_dist_best = fin_size_dist_best[k1:k2]
        #ax[1,1].plot(sizes*1e4, size_dist*(sizes**4), label=component.name)
        ax[1,1].plot(sizes*1e4, size_dist_start*(sizes**4), pfmts[i])
        ax[1,1].plot(sizes*1e4, size_dist_best*(sizes**4), pfmts_best[i])
        ax[1,1].errorbar(sizes*1e4, size_dist*(sizes**4),yerr=[fin_size_dist_munc[k1:k2]*(sizes**4), fin_size_dist_punc[k1:k2]*(sizes**4)],
                         label=component.name, fmt=pfmts_fin[i])
        k1 += component.n_sizes
    ax[1,1].legend(loc=2)

    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    ax[1,1].set_xlabel('a [$\mu$m]')
    ax[1,1].set_ylabel('$a^4 N_d(a)/N(HI)$')

    fig.savefig('dgfit_test.png')

    # save the final size distributions
    basename = 'dgfit_test'
    dustmodel.save(basename + '_sizedist.fits')
    
    if args.walkers:
        # plot the walker chains for all parameters
        fig, ax = pyplot.subplots(ndim/5, sharex=True, figsize=(13,13))
        walk_val = np.arange(nsteps)
        for i in range(1,ndim,5):
            for k in range(nwalkers):
                ax[i/5].plot(walk_val,sampler.chain[k,:,i],'-')
                #ax[i].set_ylabel(var_names[i])
        fig.savefig(basename+'_walker_param_values.png')

    if args.triangle:

        # plot the 1D and 2D likelihood functions in a traditional triangle plot
        #fig = triangle.corner(samples, labels=var_names, show_titles=True, extents=param_limits)
        fig = triangle.corner(samples)
        fig.savefig(basename+'_param_triangle.png')

