#
# File name:   make_param_pickle.py
# Date:        2019/01/10 14:36
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np
import pickle

# Density function parameters
rho_w_a = [-0.420429e-1, 0.518217702, 0.5627208e-1, 0.344164179e-1]
rho_w_a_old = [-0.420429e1, 0.518217702, 0.5627208e-1, 0.344164179e-1]

rho_w_r = [2.5, 3.1, 3.5, 4.9]

# Embedding function parameters
F_w_a = [-5.946454, -0.049477]
Fc_w_a = [-5.524855802, 2.317313103e-1, -3.665345949e-2, 8.989367404e-3]

# Pair potential parameters
arr = np.array([
(1,  0.960851701343041e2, 2.5648975),
(2, -0.184410923895214e3, 2.6297950),
(3,  0.935784079613550e2, 2.6946925),
(4, -0.798358265041677e1, 2.8663175),
(5,  0.747034092936229e1, 2.9730450),
(6, -0.152756043708453e1, 3.0797725),
(7,  0.125205932634393e1, 3.5164725),
(8,  0.163082162159425e1, 3.8464450),
(9, -0.141854775352260e1, 4.1764175),
(10,-0.819936046256149e0, 4.7008450),
(11, 0.198013514305908e1, 4.8953000),
(12,-0.696430179520267e0, 5.0897550),
(13, 0.304546909722160e-1,5.3429525),
(14,-0.163131143161660e1, 5.4016950),
(15, 0.138409896486177e1, 5.4604375)])

V_w_a = list(arr[:,1])
V_w_r = list(arr[:,2])


if __name__ == "__main__":

    # make a list of dicts of parameters with corresponding hyperparameters
    hyperparams = {}
    hyperparams['pair'] = V_w_r
    hyperparams['edens'] = rho_w_r
    params = {}
    params['pair'] = V_w_a
    params['edens'] = rho_w_a
    params['embed'] = F_w_a

    par_dict = {'hyperparams':hyperparams, 'params': params}

    par_list = [par_dict]

    with open('marinica_params.pickle', 'wb') as fo:
        pickle.dump(par_list, fo, protocol=2)
