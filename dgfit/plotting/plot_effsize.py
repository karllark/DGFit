# plot the effective grain size as a function of wavelength for a specific
# dust grain model

import argparse
import numpy as np
import matplotlib.pyplot as plt
import importlib.resources as importlib_resources

from dgfit.dustmodel import DustModel, WDDustModel

# from dgfit.dustgrains import DustGrains


def main():

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--everynth", type=int, default=1, help="Use every nth grain size"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get the dust model on the full wavelength grid
    compnames = ["astro-silicates", "astro-carbonaceous"]
    ref = importlib_resources.files("dgfit") / "data"
    with importlib_resources.as_file(ref) as data_path:
        dustmodel_full = DustModel(
            componentnames=compnames,
            path=str(data_path) + "/indiv_grain/",
            every_nth=args.everynth,
        )

    # WD model
    dustmodel = WDDustModel(dustmodel=dustmodel_full)

    plot(dustmodel, args.png, args.pdf, args.eps)


def plot(dustmodel, png=False, pdf=False, eps=False):

    # set size distributions
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
    cparams = dustmodel.parameters["Radiation field"]
    p0 += [cparams["RF"]]
    dustmodel.set_size_dist(p0)

    # compute the contribution of each grain to the total
    _effcabs_all = np.empty(
        (component.n_wavelengths, dustmodel.n_components, component.n_sizes - 1)
    )
    _effcsca_all = np.empty(
        (component.n_wavelengths, dustmodel.n_components, component.n_sizes - 1)
    )
    _effsize = np.empty((dustmodel.n_components, component.n_sizes - 1))
    _effcabs_sum = np.empty((component.n_wavelengths, dustmodel.n_components))
    _effcsca_sum = np.empty((component.n_wavelengths, dustmodel.n_components))
    _emission_sum = np.empty((component.n_wavelengths, dustmodel.n_components))

    _emission_all = np.empty(
        (component.n_wavelengths, dustmodel.n_components, component.n_sizes - 1)
    )
    for k, component in enumerate(dustmodel.components):

        # do a very simple integration (later this could be made more complex)
        deltas = 0.5 * (
            component.sizes[1 : component.n_sizes]
            - component.sizes[0 : component.n_sizes - 1]
        )
        sizedist1 = component.size_dist[0 : component.n_sizes - 1]
        sizedist2 = component.size_dist[1 : component.n_sizes]
        interpolated_emission = component.interpol_emission(1)
        for i in range(component.n_wavelengths):
            _effcabs_all[i, k, :] = deltas * (
                (component.cabs[0 : component.n_sizes - 1, i] * sizedist1)
                + (component.cabs[1 : component.n_sizes, i] * sizedist2)
            )
            _effcsca_all[i, k, :] = deltas * (
                (component.csca[0 : component.n_sizes - 1, i] * sizedist1)
                + (component.csca[1 : component.n_sizes, i] * sizedist2)
            )
            _emission_all[i, k, :] = deltas * (
                (interpolated_emission[0 : component.n_sizes - 1, i] * sizedist1)
                + (interpolated_emission[1 : component.n_sizes, i] * sizedist2)
            )

        # *not* faster to use numexpr (tested in 2015)
        _effcabs_sum[:, k] = np.sum(_effcabs_all[:, k, :], axis=1)
        _effcsca_sum[:, k] = np.sum(_effcsca_all[:, k, :], axis=1)
        _emission_sum[:, k] = np.sum(_emission_all[:, k, :], axis=1)

        _effsize[k, :] = (
            0.5
            * 1e4
            * (
                component.sizes[1 : component.n_sizes]
                + component.sizes[0 : component.n_sizes - 1]
            )
        )

    _effcext_sum = _effcabs_sum + _effcsca_sum
    _effcext_sum_allcomp = np.sum(_effcext_sum, axis=1)
    _emission_sum_allcomp = np.sum(_emission_sum, axis=1)

    # setup the plots
    fontsize = 16
    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, fax = plt.subplots(ncols=2, nrows=2, figsize=(14, 9))

    symvals = ["b--", "g--"]

    # plot size distributions
    ax = fax[1, 1]
    for k, component in enumerate(dustmodel.components):
        ax.plot(
            component.sizes * 1e4, component.size_dist, symvals[k], label=component.name
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-12, 1e2)
    ax.set_xlabel(r"a $[\mu m]$", fontsize=fontsize)
    ax.set_ylabel(r"$n_H^{-1} dn/da$", fontsize=fontsize)
    ax.set_title("Size distributions")
    ax.legend()

    # plot extinction
    ax = fax[0, 1]
    for k, component in enumerate(dustmodel.components):
        ax.plot(component.wavelengths, _effcext_sum[:, k], symvals[k])
    ax.plot(component.wavelengths, _effcext_sum_allcomp, "k-", label="total")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-25, 1e-20)
    ax.set_xlabel(r"$\lambda$ $[\mu m]$", fontsize=fontsize)
    ax.set_ylabel(r"$A(\lambda)/N(H)$", fontsize=fontsize)
    ax.set_title("Extinction")
    ax.legend()

    # plot the contribution of each grains by component to the total
    ax = fax[0, 0]
    esize = np.empty((component.n_wavelengths, dustmodel.n_components))
    for k, component in enumerate(dustmodel.components):
        for i in range(component.n_wavelengths):
            weights = (
                _effcabs_all[i, k, :] + _effcsca_all[i, k, :]
            ) / _effcext_sum_allcomp[i]
            esize[i, k] = np.sum(_effsize[k, :] * weights) / np.sum(weights)

        ax.plot(component.wavelengths, esize[:, k], symvals[k], label=component.name)

    esize_all = np.empty((component.n_wavelengths))
    for i in range(component.n_wavelengths):
        weights = (
            _effcabs_all[i, :, :] + _effcsca_all[i, :, :]
        ) / _effcext_sum_allcomp[i]
        esize_all[i] = np.sum(_effsize * weights) / np.sum(weights)
    ax.plot(component.wavelengths, esize_all, "k-", label="Total")

    ax.set_xscale("log")
    ax.set_ylim(0.0, 0.4)
    ax.set_xlabel(r"$\lambda$ $[\mu m]$", fontsize=fontsize)
    ax.set_ylabel(r"average a $[\mu m]$", fontsize=fontsize)
    ax.set_title("Average grain size for extinction")

    ax = fax[1, 0]
    esize = np.empty((component.n_wavelengths, dustmodel.n_components))
    for k, component in enumerate(dustmodel.components):
        for i in range(component.n_wavelengths):
            weights = _emission_all[i, k, :] / _emission_sum_allcomp[i]
            esize[i, k] = np.sum(_effsize[k, :] * weights) / np.sum(weights)

        ax.plot(component.wavelengths, esize[:, k], symvals[k], label=component.name)

    esize_all = np.empty((component.n_wavelengths))
    for i in range(component.n_wavelengths):
        weights = _emission_all[i, :, :] / _emission_sum_allcomp[i]
        esize_all[i] = np.sum(_effsize * weights) / np.sum(weights)
    ax.plot(component.wavelengths, esize_all, "k-", label="Total")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\lambda$ $[\mu m]$", fontsize=fontsize)
    ax.set_ylabel(r"average a $[\mu m]$", fontsize=fontsize)
    ax.set_title("Average grain size for emission")

    fig.tight_layout()

    # show or save
    basename = "effsize"
    if png:
        fig.savefig(basename + ".png")
    elif eps:
        fig.savefig(basename + ".eps")
    elif pdf:
        fig.savefig(basename + ".pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()
