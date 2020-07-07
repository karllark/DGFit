#!/usr/bin/env python

from __future__ import print_function

import sys
import time
import argparse

import numpy as np

import emcee

from dgfit.dustmodel import DustModel, MRNDustModel, WDDustModel
from dgfit.obsdata import ObsData


def DGFit_cmdparser():

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizedisttype",
        default="MRN",
        choices=["bins", "MRN", "WD"],
        help="Size distribution type",
    )
    parser.add_argument(
        "--fitobs",
        nargs="+",
        default="all",
        choices=["extinction", "iremission", "abundance", "albedo", "g", "all"],
        help="Which observations to fit",
    )
    parser.add_argument(
        "-f",
        "--fast",
        help="Use minimal walkers, steps, burns to debug code",
        action="store_true",
    )
    parser.add_argument(
        "-s", "--slow", help="Use lots of walkers, n_steps, n_burn", action="store_true"
    )
    parser.add_argument(
        "--nburn", type=int, default=500, help="Number of samples for burn"
    )
    parser.add_argument(
        "--nsteps", type=int, default=1000, help="Number of samples for full run"
    )
    parser.add_argument(
        "--everynth", type=int, default=5, help="Use every nth grain size"
    )
    parser.add_argument(
        "--chain", action="store_true", help="Store the gain in an ascii file"
    )
    parser.add_argument(
        "--limit_abund", action="store_true", help="Limit based on abundances"
    )
    parser.add_argument(
        "--usemin",
        action="store_true",
        help="Find min before EMCEE (does not work yet)",
    )
    parser.add_argument(
        "-r", "--read", default=None, help="Read size distribution from disk"
    )
    parser.add_argument(
        "-t", "--tag", default="dgfit_test", help="basename to use for output files"
    )
    parser.add_argument(
        "-c", "--cpus", metavar=int, default=4, help="number of cpus to use"
    )
    parser.add_argument(
        "--nolarge", action="store_true", help="Deweight a > 0.5 micron by 1e-10"
    )
    parser.add_argument("--smc", help="use an SMC sightline", action="store_true")

    return parser


def set_obs_for_fitting(obdata, fitobs):
    """
    parse the requested list of observations for fitting and set the
    appropriate variables
    """

    fitobs_list = []
    if not (obsdata.fit_extinction and (("extinction" in fitobs) or ("all" in fitobs))):
        obsdata.fit_extinction = False
    else:
        fitobs_list.append("extinction")
    if not (obsdata.fit_abundance and (("abundance" in fitobs) or ("all" in fitobs))):
        obsdata.fit_abundance = False
    else:
        fitobs_list.append("abundance")
    if not (
        obsdata.fit_ir_emission and (("iremission" in fitobs) or ("all" in fitobs))
    ):
        obsdata.fit_ir_emission = False
    else:
        fitobs_list.append("ir_emission")
    if not (obsdata.fit_scat_a and (("albedo" in fitobs) or ("all" in fitobs))):
        obsdata.fit_scat_a = False
    else:
        fitobs_list.append("scat albedo")
    if not (obsdata.fit_scat_g and (("g" in fitobs) or ("all" in fitobs))):
        obsdata.fit_scat_g = False
    else:
        fitobs_list.append("scat g")

    return fitobs_list


if __name__ == "__main__":

    parser = DGFit_cmdparser()

    args = parser.parse_args()

    # set the basename of the output
    basename = "%s_%s" % (args.tag, args.sizedisttype)

    # save the start time
    start_time = time.clock()

    # emcee parameters
    if args.fast:
        print("using the fast params")
        nsteps = 100
        burn = 50
    elif args.slow:
        print("using the slow params")
        nsteps = 10000
        burn = 5000
    else:
        burn = int(args.nburn)
        nsteps = int(args.nsteps)

    # get the observed data
    if args.smc:
        path = "dgfit/data/smc_azv215"
        obsdata = ObsData(
            "%s/azv215_50p_ext.fits" % path,
            "%s/azv215_avnhi.dat" % path,
            "%s/SMC_AzV215_abundances.dat" % path,
            None,
            None,
        )
    else:
        path = "dgfit/data/mw_rv31"
        obsdata = ObsData(
            [
                "%s/MW_diffuse_Gordon09_band_ext.dat" % path,
                "%s/MW_diffuse_Gordon09_iue_ext.dat" % path,
                "%s/MW_diffuse_Gordon09_fuse_ext.dat" % path,
            ],
            "%s/MW_diffuse_Gordon09_avnhi.dat" % path,
            "%s/MW_diffuse_Jenkins09_abundances.dat" % path,
            "%s/MW_diffuse_Compiegne11_ir_emission.dat" % path,
            "%s/dust_scat.dat" % path,
            ext_tags=["band", "iue", "fuse"],
            scat_path="%s/Scat_Data/" % path,
        )

    # determine what to fit based on what exists and the commandline args
    fitobs_list = set_obs_for_fitting(obsdata, args.fitobs)

    # get the dust model on the full wavelength grid
    compnames = ["astro-silicates", "astro-carbonaceous"]
    dustmodel_full = DustModel(
        componentnames=compnames,
        path="dgfit/data/indiv_grain/",
        every_nth=args.everynth,
    )

    sizedisttype = args.sizedisttype
    if sizedisttype == "MRN":
        # define the fitting model
        dustmodel = MRNDustModel(dustmodel=dustmodel_full, obsdata=obsdata)

        # initial guesses at parameters
        #    starting range is 0.001 to 1 micron
        p0 = []
        for component in dustmodel.components:
            cparams = dustmodel.parameters[component.name]
            p0 += [cparams["C"], cparams["alpha"], cparams["a_min"], cparams["a_max"]]

        # need to set dust model size distribution
        dustmodel.set_size_dist(p0)

    elif sizedisttype == "WD":
        dustmodel = WDDustModel(dustmodel=dustmodel_full, obsdata=obsdata)

        # initial guesses at parameters
        p0 = []
        for component in dustmodel.components:
            if component.name == "astro-silicates":
                cparams = dustmodel.parameters["astro-silicates"]
                p0 += [
                    cparams["C_s"],
                    cparams["a_ts"],
                    cparams["alpha_s"],
                    cparams["beta_s"],
                ]
            else:
                cparams = dustmodel.parameters["astro-carbonaceous"]
                p0 += [
                    cparams["C_g"],
                    cparams["a_tg"],
                    cparams["alpha_g"],
                    cparams["beta_g"],
                    cparams["a_cg"],
                    cparams["b_C"],
                ]

        # need to set dust model size distribution
        dustmodel.set_size_dist(p0)

    elif sizedisttype == "bins":
        dustmodel = DustModel(dustmodel=dustmodel_full, obsdata=obsdata)

        # replace the default size distribution with one from a file
        if args.read is not None:
            dustmodel.read_sizedist_from_file(args.read)

        else:
            # check that the default size distributions give approximately
            #     the right level of the A(lambda)/N(HI) curve
            # if not, adjust the overall level of the size distributions to
            #     get them close
            results = dustmodel.eff_grain_props(obsdata)
            cabs = results["cabs"]
            csca = results["csca"]
            dust_alnhi = 1.086 * (cabs + csca)
            ave_model = np.average(dust_alnhi)
            ave_data = np.average(obsdata.ext_alnhi)
            ave_ratio = ave_data / ave_model
            if (ave_ratio < 0.5) | (ave_ratio > 2):
                for component in dustmodel.components:
                    component.size_dist *= ave_ratio

                    # deweight large grains (test)
                    if args.nolarge:
                        (indxs,) = np.where(component.sizes > 0.5e-4)
                        if len(indxs) > 0:
                            print("deweighting sizes > 0.5 micron")
                            component.size_dist[indxs] *= 1e-10

        # inital guesses at parameters
        p0 = dustmodel.components[0].size_dist
        for k in range(1, dustmodel.n_components):
            p0 = np.concatenate([p0, dustmodel.components[k].size_dist])

    else:
        print("Size distribution choice not known")
        exit()

    # save the starting model
    dustmodel.save_results(basename + "_sizedist_start.fits", obsdata)

    # setup time
    setup_time = time.clock()
    print("setup time taken: ", (setup_time - start_time) / 60.0, " min")

    # more emcee setup
    ndim = len(p0)
    nwalkers = 2 * ndim

    print("fitting ", fitobs_list)
    print("# params = %i" % ndim)
    print("# walkers = %i" % nwalkers)
    print("# burn = %i" % burn)
    print("# steps = %i" % nsteps)

    # setting up the walkers to start "near" the inital guess
    p = dustmodel.initial_walkers(p0, nwalkers)

    # for pc in p:
    #    print(dgfit_model.lnprob(pc, obsdata, dustmodel))
    # exit()

    # setup the sampler
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        dustmodel.lnprob,
        args=(obsdata, dustmodel),
        threads=int(args.cpus),
    )

    # incrementally save the full chain (burn+run) to a file
    if args.chain:
        inc_prog_fname = "%s_chain.dat" % basename
        f = open(inc_prog_fname, "w")
        f.close()

    # burn in the walkers
    # pos, prob, state = sampler.run_mcmc(p, burn)
    print("burn")
    width = 60
    for i, result in enumerate(sampler.sample(p, iterations=burn)):
        n = int((width + 1) * float(i) / burn)
        sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (width - n)))

        if args.chain:
            position = result[0]
            f = open(inc_prog_fname, "a")
            for k in range(position.shape[0]):
                pos_str = " ".join(str(e) for e in position[k])
                f.write("{0:4d} {1:s}\n".format(k, pos_str))
            f.close()

    sys.stdout.write("\n")

    # decompose so the next sampling is done from the last position
    pos, prob, state = result

    # reset the sampler
    sampler.reset()

    # do the full sampling
    print("afterburn")
    width = 60
    save_frac = 0.1
    if nsteps > 1000:
        save_frac = 0.025
    targ_out = int(save_frac * nsteps)
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, rstate0=state)):
        n = int((width + 1) * float(i) / nsteps)
        sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (width - n)))

        # incrementally save the chain
        if args.chain:
            position = result[0]
            f = open(inc_prog_fname, "a")
            for k in range(position.shape[0]):
                pos_str = " ".join(str(e) for e in position[k])
                f.write("{0:4d} {1:s}\n".format(k, pos_str))
            f.close()

        # output the size distribution
        if i > targ_out:
            oname = "%s_sizedist_%i.fits" % (basename, targ_out)
            dustmodel.save_50percentile_results(oname, sampler, obsdata, cur_step=i)
            oname = "%s_sizedist_best_%i.fits" % (basename, targ_out)
            dustmodel.save_best_results(oname, sampler, obsdata, cur_step=i)
            targ_out += int(save_frac * nsteps)

    sys.stdout.write("\n")

    emcee_time = time.clock()
    print("emcee time taken: ", (emcee_time - setup_time) / 60.0, " min")

    # best fit dust params
    oname = "%s_sizedist_best_fin.fits" % (basename)
    dustmodel.save_best_results(oname, sampler, obsdata)

    # 50p dust params
    oname = "%s_sizedist_fin.fits" % (basename)
    dustmodel.save_50percentile_results(oname, sampler, obsdata)
