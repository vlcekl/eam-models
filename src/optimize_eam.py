from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division
try:
    xrange = xrange
    # We have Python 2
except:
    xrange = range
    # We have Python 3


import os
import sys
import re
import numpy as np
import pickle
from scipy.optimize import fmin
import multiprocessing as mp
import functools


# statmech library setup
sys.path.append('../../statmechlib')
from statmechlib.preprocessing import Trajectory, select_nodes, scale_configuration, pair_dist_cutoff
from statmechlib.forcefields import sd2_loss, utot_EAM_per_atom, utot_EAM_per_box, ftot_EAM, udif_print, u_core
from statmechlib.preprocessing import universal_eos
from statmechlib.preprocessing import pair_dist, force_targ, get_stats_EAM_per_atom, get_stats_EAM_per_box

target_raw = '../data/target_raw'
target_proc = '../data/target_processed'
working = '../data/working'

#stats_file = 'stats_marinica' # 'stats_samples'
#target_file = 'target_marinica' # 'target_samples'
stats_file = 'stats_manyknots' # 'stats_samples'
target_file = 'target_manyknots' # 'target_samples'
params_file = 'marinica_params'


with open(os.path.join(working, stats_file+'.pickle'), 'rb') as fi:
    stats_all = pickle.load(fi)

with open(os.path.join(working, target_file+'.pickle'), 'rb') as fi:
    targets = pickle.load(fi)

with open(os.path.join(working, params_file + '.pickle'), 'rb') as fi:
    param_list = pickle.load(fi)  


def find_index(select_list, full_list):
    knots = []
    for sel in select_list:
        for i, elem in enumerate(full_list):
            if abs(sel - elem) < 1e-9:
                knots.append(i)
                break
    
    assert len(knots) == len(select_list), "Knots and select_list lengths do not match"
    
    return knots


# find indices for chosen spline knots (this one is for Marinica potential)
#pair_knots = [ 2.5648975,  2.629795 ,  2.6946925,  2.8663175,  2.973045 ,
#        3.0797725,  3.5164725,  3.846445 ,  4.1764175,  4.700845 ,
#        4.8953   ,  5.089755 ,  5.3429525,  5.401695 ,  5.4604375]
#edens_knots = [ 2.5,  3.1,  3.5,  4.9]


# extended knots
pair_knots = [ 2.45, 2.5648975,  2.629795 ,  2.6946925,  2.8663175,  2.973045 ,
        3.0797725,  3.5164725,  3.846445 ,  4.1764175,  4.700845 ,
        4.8953   ,  5.089755 ,  5.3429525,  5.401695 ,  5.4604375, 5.78]
edens_knots = [ 2.5,  3.1,  3.5,  4.9, 5.5]

pair_index = find_index(pair_knots, stats_all['hyperparams']['pair'])
edens_index = find_index(edens_knots, stats_all['hyperparams']['edens'])
print('pair_index:', pair_index)
print('edens_index:', edens_index)


# select spline knots
#pair_index = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 15, 16, 17, 18]
#edens_index = [0, 7, 8, 14]

p_ix = np.array([True if i in pair_index else False for i in range(len(stats_all['hyperparams']['pair']))])
m_ix = np.array([True if i in edens_index else False for i in range(len(stats_all['hyperparams']['edens']))])

stats = select_nodes(stats_all, p_ix, m_ix)

print("pair:", np.array(stats_all['hyperparams']['pair'])[p_ix])
print("edens:", np.array(stats_all['hyperparams']['edens'])[m_ix])
print('pars', param_list[0]['hyperparams'])


multi_pars = []
for params in param_list:
    eam_params = list(params['params']['embed']) + [0.0] = list(params['params']['pair']) + [0.0] + list(params['params']['edens']) + [0.0]
    multi_pars.append(np.array(eam_params))


multi_pars.append(np.array([ -6.67892458e+00,  -6.69057233e-02, 0.0, -7.88558288e+02,
         8.64431119e+02,  -3.22117669e+02,   1.66223488e+01,
         1.61253439e+01,  -7.08057431e+00,   3.36500438e-01,
         1.23818786e+00,  -9.77125576e-01,  -7.10990440e-01,
         1.96923157e+00,  -7.11735164e-01,   3.61922715e-02,
        -1.61997750e+00,   1.35041921e+00,  0.0, -3.71430662e+01,
         8.77780820e-01,  -5.71131123e-02,   3.98140807e-02, 0.0]))

multi_pars.append(np.array([ -9.42391613e+00,  -4.24391684e-02, 0.0, -2.94587357e+02,
         2.70972700e+02,  -1.04787033e+02,  -2.98238182e+00,
        -4.35229417e-01,   2.92307252e+01,  -3.93677905e+00,
         3.09415379e+00,  -2.40064822e+00,  -8.22415578e-01,
         1.49304711e+00,  -3.38081838e-02,   3.54435455e-01,
        -1.76673992e+00,   1.15161862e+00,  0.0, -3.15350855e+01,
         4.07693197e+00,  -7.52508632e-01,   3.10740894e-02, 0.0]))

multi_pars.append(np.array([ -8.15036857e+00,  -1.34492508e-01,  0.0, -6.40562492e+02,
         5.60045606e+02,  -1.62774419e+02,  -2.42163837e+01,
         4.61349632e+01,  -1.21579989e+01,   6.09852042e-01,
         1.39564841e+00,  -1.26608655e+00,  -7.44806685e-01,
         2.00404873e+00,  -6.91376693e-01,   2.17652121e-02,
        -1.61632910e+00,   1.37087609e+00,  0.0, -3.76175083e+01,
         1.56626002e+00,  -1.74622706e-01,   3.14905374e-02, 0.0]))

multi_pars.append(np.array([ -5.10420932e+00,  -3.58550988e-03, 7.40616359e-01,
        -2.51494375e+02,   2.13567082e+02,  -5.65076897e+01,
        -3.54179888e+01,   2.09710220e+01,   2.12912000e+01,
        -3.85412517e+00,   3.97807318e+00,  -2.40631602e+00,
        -9.66381896e-01,   1.51743634e+00,  -1.46796026e-01,
         3.69930757e-01,  -1.74239343e+00,   1.16099701e+00,
        -4.62036585e-03,  -5.70390094e+01,   1.02379444e+01,
        -1.57673128e+00,   3.13141586e-02,   1.28354461e-02]))


#targets['relax']['energy'], len(stats['relax']['energy']), len(targets['relax']['temp'])#, stats['relax']['energy'][-2:]
assert len(targets['relax']['temp']) == 7, "Temperature wrong"
assert len(stats['relax']['energy']) == 7, "Stats energy wrong"
#targets['relax']['temp'] = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
#stats['relax']['energy'] = stats['relax']['energy'][:-1]
#stats['relax']['forces'] = stats['relax']['forces'][:-1]

#multi_pars.append(np.array([-6.25140612e-01,  5.84924686e-05,  1.50463356e+01, -2.45628080e+00,
#        2.75467959e+00, -1.82521768e+00,  5.82623873e-01, -8.10349686e-02, 0.0]))
#for _ in range(n_walkers):
#    pars += np.random.standard_normal(multi_pars[0].shape)*0.01


weight_dict = {'md':10.0, 'relax':10.0, 'eos_bcc':1.0, 'eos_fcc':0.4}

for key in targets:
    targets[key]['weight'] = weight_dict[key]


def optimize_EAM_mp(pars, targets, stats, utot_func):
    
    output = fmin(sd2_loss, pars, args=(targets, stats, utot_EAM_per_atom, None, 0.05, 1), maxiter=100000, maxfun=100000, disp=0, full_output=1,ftol=1e-6)

    return tuple([output[1], output[0]])


params_output = []

pool = mp.Pool()

for ieam in range(1):#5, len(stats_opts[it]['hyperparams'])):

    get_sd = functools.partial(sd2_loss, targets=targets, stats=stats, utot_func=utot_EAM_per_atom, ftot_func=None, dl=0.05, verbose=0)

    get_opt = functools.partial(optimize_EAM_mp, targets=targets, stats=stats, utot_func=utot_EAM_per_atom)

    # initial ordering
    print('Initial params:', multi_pars)
    
    distances = pool.map(get_sd, multi_pars)
    
    optimal_params = zip(distances, multi_pars)
    best_params = sorted(optimal_params, key=lambda param: param[0], reverse=True)
    m_pars = [p[1] for p in best_params]

    print('Best params:', best_params)
    print('ieam:', ieam)
    #m_pars = multi_pars

    for i in range(8):
        optimal_params = pool.map(get_opt, m_pars)
        best_params = sorted(optimal_params, key=lambda param: param[0], reverse=True)
        m_pars = [p[1] for p in best_params]
        
        print("Iteration {}, best params: {}".format(i, best_params))
        sys.stdout.flush()

    print('Final Best params:', best_params)
    
    params_output.append(best_params)

pool.close()
pool.join()


