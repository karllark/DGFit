import pkg_resources
import time
import argparse

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import minimize
import emcee
from multiprocessing import Pool
import corner

from dgfit.dustmodel import DustModel, MRNDustModel, WDDustModel
from dgfit.obsdata import ObsData


def DGFit_cmdparser():
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizedisttype",
        default="WD",
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
        "--mcmc", help="Do MCMC sampling using emcee package", action="store_true"
    )

    parser.add_argument(
        "-f",
        "--fast",
        help="MCMC: Use minimal walkers, steps, burns to debug code",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--slow",
        help="MCMC: Use lots of walkers, n_steps, n_burn",
        action="store_true",
    )
    parser.add_argument(
        "--burnfrac",
        type=int,
        default=500,
        help="Fractinal portion of nsteps for burn in",
    )
    parser.add_argument(
        "--nsteps", type=int, default=1000, help="Number of samples for full run"
    )
    parser.add_argument(
        "--everynth", type=int, default=2, help="Use every nth grain size"
    )
    parser.add_argument(
        "--chain", action="store_true", help="Store the chain in an ascii file"
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

    return parser


def set_obs_for_fitting(obsdata, fitobs):
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


def main():
    parser = DGFit_cmdparser()

    args = parser.parse_args()

    # set the basename of the output
    basename = f"{args.tag}_{args.sizedisttype}"

    # save the start time
    start_time = time.process_time()

    # emcee parameters
    if args.fast:
        print("using the fast params")
        nsteps = 100
        burnfrac = 0.1
    elif args.slow:
        print("using the slow params")
        nsteps = 10000
        burnfrac = 0.2
    else:
        burnfrac = int(args.burnfrac)
        nsteps = int(args.nsteps)

    # get the location of the provided data
    data_path = pkg_resources.resource_filename("dgfit", "data/")

    # get the observed data
    path = f"{data_path}/mw_rv31"
    obsdata = ObsData(
        #        [
        #            f"{path}/MW_diffuse_Gordon09_band_ext.dat",
        #            f"{path}/MW_diffuse_Gordon09_iue_ext.dat",
        #            f"{path}/MW_diffuse_Gordon09_fuse_ext.dat",
        #        ],
        [f"{path}/MW_diffuse_Gordon23_ext.dat"],
        f"{path}/MW_diffuse_Gordon09_avnhi.dat",
        f"{path}/MW_diffuse_Jenkins09_abundances.dat",
        f"{path}/MW_diffuse_Compiegne11_ir_emission.dat",
        f"{path}/dust_scat.dat",
        ext_tags=["band", "iue", "fuse"],
        scat_path=f"{path}/Scat_Data/",
    )

    # determine what to fit based on what exists and the commandline args
    fitobs_list = set_obs_for_fitting(obsdata, args.fitobs)

    # get the dust model on the full wavelength grid
    compnames = ["astro-silicates", "astro-carbonaceous"]
    dustmodel_full = DustModel(
        componentnames=compnames,
        path=f"{data_path}/indiv_grain/",
        every_nth=args.everynth,
    )

    print(f"# of grain sizes = {len(dustmodel_full.components[0].sizes)}")

    sizedisttype = args.sizedisttype
    pnames = []
    if sizedisttype == "MRN":
        # define the fitting model
        dustmodel = MRNDustModel(dustmodel=dustmodel_full, obsdata=obsdata)

        # initial guesses at parameters
        #    starting range is 0.001 to 1 micron
        p0 = []
        for component in dustmodel.components:
            cparams = dustmodel.parameters[component.name]
            p0 += [cparams["C"], cparams["alpha"], cparams["a_min"], cparams["a_max"]]
            pnames += cparams.keys()

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
            pnames += cparams.keys()

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

        # initial guesses at parameters
        p0 = []
        for k in range(0, dustmodel.n_components):
            p0 = np.concatenate([p0, dustmodel.components[k].size_dist])
            pnames += [
                f"c{k+1}_s{kk}" for kk in range(len(dustmodel.components[k].size_dist))
            ]

    else:
        print("Size distribution choice not known")
        exit()

    # save the starting model
    dustmodel.save_results(basename + "_sizedist_start.fits", obsdata)

    # setup time
    setup_time = time.process_time()
    print("setup time taken: ", (setup_time - start_time) / 60.0, " min")

    # do simple optimization to find the best fit
    def nll(*args):
        return -dustmodel.lnprob(*args)

    soln = minimize(nll, p0, args=(obsdata, dustmodel), method="Nelder-Mead")
    opt_params = soln.x
    dustmodel.set_size_dist_parameters(opt_params)

    oname = f"{basename}_sizedist_best_optimizer.fits"
    # TODO: add saving of the size distribution parameters for the analytic forms
    dustmodel.save_results(oname, obsdata)

    opt_time = time.process_time()
    print("optimizer time taken: ", (opt_time - setup_time) / 60.0, " min")

    if args.mcmc:
        p0 = opt_params
        # more emcee setup
        ndim = len(p0)
        nwalkers = 2 * ndim

        print(f"fitting {fitobs_list}")
        print(f"# params = {ndim}")
        print(f"# walkers = {nwalkers}")
        print(f"# burnfrac = {burnfrac}")
        print(f"# steps = {nsteps}")

        # setting up the walkers to start "near" the inital guess
        p = dustmodel.initial_walkers(p0, nwalkers)

        # for pc in p:
        #    print(dgfit_model.lnprob(pc, obsdata, dustmodel))
        # exit()

        # Set up the backend to save the samples for the emcee runs
        emcee_samples_file = f"{basename}_chain.h5"
        backend = emcee.backends.HDFBackend(emcee_samples_file)
        backend.reset(nwalkers, ndim)

        # setup the sampler
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(
                nwalkers,
                ndim,
                dustmodel.lnprob,
                args=(obsdata, dustmodel),
                pool=pool,
                backend=backend,
            )

            # do the sampling
            sampler.run_mcmc(p, nsteps, progress=True)

        emcee_time = time.process_time()
        print("emcee time taken: ", (emcee_time - opt_time) / 60.0, " min")

        # plot the walker chains for all parameters
        nwalkers, nsteps, ndim = sampler.chain.shape
        fig, ax = plt.subplots(ndim, sharex=True, figsize=(13, 13))
        walk_val = np.arange(nsteps)
        for i in range(ndim):
            for k in range(nwalkers):
                ax[i].plot(walk_val, sampler.chain[k, :, i], "-")
                ax[i].set_ylabel(pnames[i])
        fig.savefig(f"{basename}_walker_param_values.png")
        plt.close(fig)

        # plot the 1D and 2D likelihood functions in a traditional triangle plot
        nwalkers, nsteps = sampler.lnprobability.shape
        # discard the 1st burn_frac (burn in)
        flat_samples = sampler.get_chain(discard=int(burnfrac * nsteps), flat=True)
        nflatsteps, ndim = flat_samples.shape
        fig = corner.corner(
            flat_samples,
            labels=pnames,
            show_titles=True,
            title_fmt=".3f",
            use_math_text=True,
        )
        fig.savefig(f"{basename}_param_triangle.png")
        plt.close(fig)

        # best fit dust params
        oname = "%s_sizedist_best_fin.fits" % (basename)
        dustmodel.save_best_results(oname, sampler, obsdata)

        # 50p dust params
        oname = "%s_sizedist_fin.fits" % (basename)
        dustmodel.save_50percentile_results(
            oname, sampler, obsdata, nburn=int(burnfrac * nsteps)
        )


if __name__ == "__main__":
    main()
