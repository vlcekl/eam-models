#!/usr/bin/python

from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division
from __future__ import unicode_literals

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
from statmechlib.preprocessing import Trajectory
from statmechlib.preprocessing import force_targ, get_stats_EAM_per_atom

#sc = list(np.round(np.linspace(2.02, 5.82, 96), 2))
sc = list(np.round(np.linspace(1.0, 5.85, 98), 2))
sc.extend([ 2.45, 2.5648975, 2.6297950, 2.6946925,
            2.8663175, 2.9730450, 3.0797725, 3.5164725,
            3.8464450, 4.1764175, 4.7008450, 4.8953000,
            5.0897550, 5.3429525, 5.4016950, 5.4604375])

sc.sort()

print(sc)

#name = 'marinica'
name_in = 'all_samples'

# Define locations of relevant datasets

working = '../data/working'
trjfile = 'trj_' + name_in + '.pickle'

name = 'pb_samples'


# read dict with trajectories information
with open(os.path.join(working, trjfile), 'rb') as fi:
    trj_fit = pickle.load(fi)


# Prepare target_data dict

#weights = {k:1.0 for k in trj_fit}
weights = {}
for key in trj_fit.keys():
    weights[key] = 1.0

target_data = {}

for key, trj in trj_fit.items():
    
    print('dataset #', key)
    sys.stdout.flush()

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
with open(os.path.join(working, "target_" + name + ".pickle"), 'wb') as fo:
    pickle.dump(target_data, fo, protocol=2)


# Prepare stats_data dict using multiprocessing

stats_data = {}

stats_data['function'] = 'EAM-cubic-spline'
stats_data['hyperparams'] = {'pair':sc, 'edens':sc}

get_stats = functools.partial(get_stats_EAM_per_atom, sc=sc, rcut=None, atom_type=None)

pool = mp.Pool()

for key, trj in trj_fit.items():
    
    print('dataset #', key)
    sys.stdout.flush()

    configs = zip(trj['xyz'], trj['box'])

    output_stats = pool.map(get_stats, configs)

    # statistics data
    stats_dict = {'energy':[], 'forces':[]}
    
    for (a1, ar, a2, ax, f1, fr, f2, c1) in output_stats:
        stats_dict['energy'].append([ar, a2, a1, ax, c1])
        stats_dict['forces'].append([fr, f2, f1])

    stats_data[key] = stats_dict
    
pool.close()
pool.join()

with open(os.path.join(working, "stats_" + name + ".pickle"), 'wb') as fo:
    pickle.dump(stats_data, fo, protocol=2)

