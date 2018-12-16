#
# # from scipy.interpolate import interp1d
# from scipy.optimize import minimize
#
# # from lmfit import minimize, Parameters
#
#     # call scipy.optimize to get a better initial guess
#     if args.usemin:
#         print(p0)
#         print(lnprob_discrete(p0, obsdata, dustmodel))
#
#         # generate the bounds
#         p0_bounds = []
#         for k in range(len(p0)):
#             p0_bounds.append((0.0, 1e20))
#
#             # setup what can be fit
#         obsdata.fit_extinction = True
#         obsdata.fit_abundance = False
#         obsdata.fit_ir_emission = False
#         obsdata.fit_scat_a = False
#         obsdata.fit_scat_g = False
#
#         # neg_lnprobsed = lambda *args: -1.0*lnprob_discrete(*args)
#         def neg_lnprobsed(*args): -1.0*lnprob_discrete(*args)
#         better_start = minimize(neg_lnprobsed, p0, args=(obsdata, dustmodel),
#                                 bounds=p0_bounds, method='L-BFGS-B')
#         print(better_start.success)
#         print(better_start.x)
#         exit()
#
#
#
#     # import scipy.optimize as op
#     # nll = lambda *args: -lnprobsed(*args)
#     # result = op.minimize(nll, p0, args=(obsdata, dustmodel))
#     # print(result)
#     # exit()
#
#     # trying with lmfit (not tested)
#     # params = Parameters()
#     # for i in range(n_params):
#     #   params.add('p'+str(i),value=p0[i],min=0.0)
#     # out = minimize(kext_residuals, params, args=(1.0/xdata, ydata))
#     # print(out.params)
#
#     # walker start check (bins model)
#         # make sure each walker starts with allowed abundances
#         if args.limit_abund:
#             for pc in p:
#                 dustmodel.set_size_dist(pc)
#                 results = dustmodel.eff_grain_props(obsdata)
#                 cabs = results['cabs']
#                 csca = results['csca']
#                 natoms = results['natoms']
#                 max_violation = 0.0
#                 for atomname in natoms.keys():
#                     cur_violation = (natoms[atomname]
#                                      / (obsdata.abundance[atomname][0] +
#                                         obsdata.abundance[atomname][1]))
#                     if cur_violation > max_violation:
#                         max_violation = cur_violation
#                 if max_violation > 2:
#                     pc *= 1.9/max_violation
#
# # could be useful....from lnprob all
#                     # hard limit at 1.5x the total possible abundaces
#                     #      (all atoms in dust)
#                     # if natoms[atomname] > 1.5*obsdata.total_abundance[atomname][0]:
#                         # print('boundary issue')
#                         # return -np.inf
#                         # pass
#                     # only add if natoms > depletions
#                     # elif natoms[atomname] > obsdata.abundance[atomname][0]:
