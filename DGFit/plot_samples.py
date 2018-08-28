from __future__ import print_function

import argparse

import numpy as np
# import matplotlib.pyplot as pyplot
# import matplotlib

import corner

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",
                        help=("file with EMCEE sampler chain"))
    args = parser.parse_args()

    samples_data = np.loadtxt(args.filename)

    nparam = len(samples_data[0, :]) - 1
    samples = samples_data[:, 1:nparam+1]

    samples = np.log10(samples[:, ::5])
    print(samples.shape)

    fig = corner.corner(samples)
    fig.savefig("%s.png" % args.filename)
