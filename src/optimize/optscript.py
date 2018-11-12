#!//anaconda/envs/py36/bin/python
#
# File name:   optscript.py
# Date:        2018/11/08 13:16
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import os
import re
import pickle
import numpy as np
from scipy.optimize import fmin

sys.path.append('../../../statmechlib')
from statmechlib.read_write import read_vasp
from statmechlib.preprocessing import Trajectory
from statmechlib.forcefields import sd2_loss, utot_EAM

if __name__ == "__main__":

    working = '../../data/working'

    # load target data
    with open(os.path.join(working, 'target_all'+'.pickle'), 'rb') as fi:
        targ_dict = pickle.load(fi)

    # load stats data
    with open(os.path.join(working, 'stats_all'+'.pickle'), 'rb') as fi:
        stats_dict = pickle.load(fi)

    # load parameters
    with open(os.path.join(working, 'pars_in'+'.pickle'), 'rb') as fi:
        pars_dict = pickle.load(fi)

    # scale energies (subtract energy of an isolated atom)
    # determine interaction energy
    print(targ_dict['relax']['energy'][0], len(targ_dict['relax']['xyz'][0]))
    u_t = targ_dict['relax']['energy'][0]/len(targ_dict['relax']['xyz'][0])
    u_e = -8.9 # external energy from atom (experimental)
    u_i = u_t - u_e # internal energy per atom (to be subtracted from all atoms)
    print(u_t, u_e, u_i)

    emin = 0.0
    esum = 0.0
    isum = 0.0
    for key, trj in targ_dict.items():
        for i in range(len(targ_dict[key]['energy'])):
            #print(key,len(targ_dict[key]['energy']), targ_dict[key]['energy'][i], targ_dict[key]['xyz'][i].shape[0])
            targ_dict[key]['energy'][i] -= u_i*targ_dict[key]['xyz'][i].shape[0]
            #enes[i] -= u_i*xyzs[i].shape[0]
            #print(key,len(targ_dict[key]['energy']), targ_dict[key]['energy'][i], targ_dict[key]['xyz'][i].shape[0])
            #print(key, targ_dict[key]['energy'][i]/targ_dict[key]['xyz'][i].shape[0])
            enex = targ_dict[key]['energy'][i]/targ_dict[key]['xyz'][i].shape[0]
            if 'fcc' in key:
                esum += enex
                isum += 1.0
                if enex < emin:
                    print(enex)
                    emin = enex

    print('Emin', emin, esum/isum)

    # temporary fix - set fcc energy at 0K to minimum energy of 300K trajectory
    targ_dict['relax']['energy'][1] = emin*targ_dict['relax']['xyz'][1].shape[0]
    targ_dict['relax']['energy'][-1] = 0.0

    print('OK energies:', targ_dict['relax']['energy'])
    targ_dict['relax']['weight'] = 50.0
    print([targ_dict[k]['weight'] for k in targ_dict])

    # prepare data for fitting
    stats = []
    target = []
    for key in list(stats_dict.keys()):
        stats.append(stats_dict[key])
        target.append(targ_dict[key])

    # prepare parameters for fitting
    pars_in = [pars_dict['embed'][0], pars_dict['embed'][2], *pars_dict['pair']]


    # fit 
#    output = fmin(sd2_loss, pars_in, args=(stats, target, utot_EAM), maxiter=100000, maxfun=100000, disp=0, full_output=1,ftol=1e-6)
#    params_uopt = output[0]
#    print("Optimized parameters:")
#    print(*params_uopt)

    multi_pars = [np.array(pars_in)]
    multi_pars.append(np.array([-1.17194534819, 4.9212636569e-05, 0.0960596087037, 16.9530837862, -1.26438173901, 2.1048867031, -0.912012405654, 0.113324291952]))

    for _ in range(10):
        pars = np.array(pars_in)
        pars += np.random.standard_normal(pars.shape)
        print(pars)
        multi_pars.append(pars)

    optimal_parameters = []
    for i, pars in enumerate(multi_pars):
        sd_ini = sd2_loss(list(pars), stats, target, utot_EAM)
        print('Initial sd2:', sd_ini)
        output = fmin(sd2_loss, list(pars), args=(stats, target, utot_EAM), maxiter=100000, maxfun=100000, disp=0, full_output=1,ftol=1e-6)
        params_uopt = output[0]
        print('Opt #', i, output[1:])
        print("Optimized parameters:")
        print(*params_uopt)
        optimal_parameters.append(tuple([output[1:], params_uopt]))

    with open(os.path.join(working, 'output_pars.pickle'), 'wb') as fo:
        pickle.dump(optimal_parameters, fo)

