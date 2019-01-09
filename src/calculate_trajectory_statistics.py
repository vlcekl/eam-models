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
import numpy as np
import pickle
import functools
import multiprocessing as mp


# statmech library setup
sys.path.append('../../statmechlib')
from statmechlib.read_write import read_vasp
from statmechlib.preprocessing import Trajectory, select_nodes, scale_configuration
from statmechlib.preprocessing import pair_dist, force_targ, get_stats_EAM_per_atom, get_stats_EAM_per_box


# Define locations of relevant datasets

target_raw = '../data/target_raw'
target_proc = '../data/target_processed'
working = '../data/working'


# read dict with trajectories information
with open(os.path.join(working, "trj_fit.pickle"), 'rb') as fi:
    trj_fit = pickle.load(fi)


# Prepare target_data dict

weights = {k:1.0 for k in trj_fit}
weights['relax'] = 10.0
weights['eos_bcc'] = 10
weights['eos_fcc'] = 0.5

target_data = {}

for key, trj in trj_fit.items():
    
    print('dataset #', key)

    # target data
    target_dict = {'type':'trajectory'}
    target_dict['weight']= weights[key]
    target_dict['box'] = trj['box']
    target_dict['xyz'] = trj['xyz']
    target_dict['energy'] = trj['energy']
    target_dict['forces'] = force_targ(trj['forces'])
    target_dict['temp'] = trj['temp']

    # save inverse temperature data (if T=0, set beta=1/300)
    target_dict['beta'] = np.empty_like(target_dict['temp'])
    for i, temp in enumerate(target_dict['temp']):
        if temp == 0.0:
            target_dict['beta'][i] = 1.0/300.0
        else:
            target_dict['beta'][i] = 1.0/temp

    target_dict['beta'] = list(target_dict['beta'])
            
    target_data[key] = target_dict

# save target data
with open(os.path.join(working, "target_samples.pickle"), 'wb') as fo:
    pickle.dump(target_data, fo)


#sc = [2.4, 2.45, 2.5, 2.65, 2.70, 2.72, 2.73, 2.74, 2.75, 2.77, 2.80, 2.85, 2.90,
#      3.0, 3.1, 3.2, 3.3, 3.45, 3.6, 3.75,
#      4.0, 4.1, 4.25, 4.356, 4.5, 4.65, 4.8,
#      5.0, 5.15, 5.3, 5.45, 5.6, 5.75]
#sc = [2.74, 3.252, 3.804, 4.356, 4.908, 5.46]

sc = list(np.linspace(2.38, 5.78, 86))

# Prepare stats_data dict using multiprocessing

stats_data = {}

stats_data['function'] = 'EAM-cubic-spline'
stats_data['hyperparams'] = {'pair':sc, 'edens':sc}

get_stats = functools.partial(get_stats_EAM_per_atom, sc=sc, rcut=None, atom_type=None)

pool = mp.Pool()

for key, trj in trj_fit.items():
    
    print('dataset #', key)

    configs = zip(trj['xyz'], trj['box'])

    output_stats = pool.map(get_stats, configs)

    # statistics data
    stats_dict = {'energy':[], 'forces':[]}
    
    for (a1, ar, a2, ax, f1, fr, f2) in output_stats:
        stats_dict['energy'].append([ar, a2, a1, ax])
        stats_dict['forces'].append([fr, f2, f1])

    stats_data[key] = stats_dict
    
pool.close()
pool.join()

with open(os.path.join(working, "stats_samples.pickle"), 'wb') as fo:
    pickle.dump(stats_data, fo)

