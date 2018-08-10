#!//anaconda/envs/py36/bin/python
#
# File name:   stats.py
# Date:        2018/07/24 16:17
# Author:      Lukas Vlcek
#
# Description: Read in configurations, calculate statistics for subsequent use
# in IFF optimization. Store data in text files
#

import sys
import os
import re
import numpy as np
import h5py
import pickle
from itertools import product
import matplotlib.pyplot as plt

# pair distances
def pair_dist(xyz, box):
    """
    Calculates nearest image pair distances between all atoms in xyz array.
    Parameters
    -----------
    xyz : numpy array
          particle x, y, z coordinates
    box : scalar or numpy array
          simulation box dimensions/shape
    Returns
    -------
    rr  : (natom, natom) numpy array of pair distances
    rx  : (natom, natom, 3) numpy array of pair distance coordinates
    """

    n_atom = xyz.shape[0] # number of atoms in a configuration
    rr = np.empty((n_atom, n_atom), dtype=float)
    rx = np.empty((n_atom, n_atom, 3), dtype=float)

    for i, pa in enumerate(xyz):
        for j, pb in enumerate(xyz):
            dp = pa - pb
            dp = np.where(dp < -0.5*box, dp + box, dp)
            dp = np.where(dp >  0.5*box, dp - box, dp)
            rr[i,j] = np.sum(dp*dp)**0.5
            rx[i,j] = dp
        
    return rr, rx

def pair_dist_minibox(xyz, box):
    """
    Calculates nearest image pair distances between all atoms in xyz array.
    Parameters
    -----------
    xyz : numpy array
          particle x, y, z coordinates
    box : scalar or numpy array
          simulation box dimensions/shape
    Returns
    -------
    rr  : (natom, natom) numpy array of pair distances
    rx  : (natom, natom, 3) numpy array of pair distance coordinates
    """

    # create box replicas in all directions
    for ib1, ib2, ib3 in product([-1, 0, 1], repeat=3):
        print('ib', ib1, ib2, ib3)

    n_atom = xyz.shape[0] # number of atoms in a configuration
    rr = np.empty((n_atom, n_atom), dtype=float)
    rx = np.empty((n_atom, n_atom, 3), dtype=float)

    for i, pa in enumerate(xyz):
        for j, pb in enumerate(xyz):
            dp = pa - pb
            dp = np.where(dp < -0.5*box, dp + box, dp)
            dp = np.where(dp >  0.5*box, dp - box, dp)
            rr[i,j] = np.sum(dp*dp)**0.5
            rx[i,j] = dp
        
    #rr = np.array(rr)
    #rx = np.array(rx)
    return rr, rx

# sufficient statistics for EAM
def get_stats_EAM(rr, rx, sc):
    """
    Takes atom pair distances and calculates sufficeint statistics needed
    for the parameterization of a cubic spline-based EAM model by Bonny et al. (2017).
 
    Parameters
    ----------
    rr : numpy array
         set of pair distances
    rx : numpy array
         set of pair distance coordinates
    sc : python list
         spline nodes

    Returns
    -------
    ar, a1, a2 : numpy arrays (len(sc))
                 atom energy-related statistics
                 el_density**0.5, el_density, el_density**2
    br, b1, b2 : numpy arrays (len(sc), natoms, 3 coordinates)
                 atom force-related statistics (gradients of energy)
                 grad(el_density**0.5), grad(el_density), grad(el_density**2)
    """
 
    n_atom = rr.shape[0]
    
    # energy-related statistics
    aa = np.empty((n_atom), dtype=float)
    ar = np.zeros((len(sc)), dtype=float)
    a1 = np.zeros_like(ar)
    a2 = np.zeros_like(ar)
    
    # force-related statistics
    br = np.zeros((len(sc), n_atom, 3), dtype=float)
    b1 = np.zeros_like(br)
    b2 = np.zeros_like(br)
    zero3 = np.zeros((3))

    # cycle over spline nodes
    for ks, rc in enumerate(sc):

        # cycle over atoms
        for i in range(n_atom):

            # sum electronic density over all neighbors of i within rc
            aa[i] = sum([(rc - r)**3 for r in rr[i] if (r < rc and r > 0.01)])

            # if el. density larger than zero, calculate force statistics
            if aa[i] > 0.0:

                # precompute a list of recurring values for force statistics
                ff = [1.5*(rc - r)**2*x/r if (r > 0.01 and r < rc) else zero3 for r, x in zip(rr[i], rx[i])]

                # sum contributions to force statistics from all neighbors of i
                b1[ks, i] = sum([2*f       for f in ff])
                br[ks, i] = sum([ -f/np.sqrt(aa[i]) for f in ff])
                b2[ks, i] = sum([4*f*aa[i] for f in ff])

        # sum contributions to energy statistics for a given spline node
        ar[ks] = np.sum(np.sqrt(aa))
        a1[ks] = np.sum(aa)
        a2[ks] = np.sum(aa**2)

    return a1, ar, a2, b1, br, b2

def force_targ(forces):

    force_flat = []
    for frc in forces:
        fr = np.concatenate((np.array([0.0]), frc.flatten(), -frc.flatten()))
        force_flat.append(fr)

    return np.array(force_flat)


if __name__ == "__main__":

    # processed data directory
    target_proc = '../../../data/target_processed'
    working = '../../../data/working'

    # Static parameters of the EAM potential: set of spline nodes
    sc = [2.56, 2.73, 3.252, 3.804, 4.20, 4.77]

    files = ['structs_0k', 'liq_4000k', 'bcc_300k']#, 'german_dft.h5']
    weights = [10.0, 1.0, 1.0]
    
    # Dictionary of statistics and target data to be used in optimization
    stats_data = {}
    target_data = {}

    # cycle over trajectory data and save target and statistics information
    for di, (filin, weight) in enumerate(zip(files, weights)):

        # Create a target dataset directory with exhaustive target information
        target_dict = {'type':'trajectory', 'weight':weight}
        with open(os.path.join(target_proc, filin+'.pickle'), 'rb') as fi:
            # read pickled trajectory dictionary
            traj_dict = pickle.load(fi)

        # save trajectory data
        target_dict['box'] = traj_dict['box']
        target_dict['xyz'] = traj_dict['xyz']
        target_dict['energy'] = traj_dict['energy']
        target_dict['temp'] = traj_dict['temp']

        # read and transform forces into (6N+1) arrays
        if 'forces' in traj_dict.keys():
            target_dict['forces'] = force_targ(traj_dict['forces'])

        # save inverse temperature data (if T=0, set beta=1/300)
        target_dict['beta'] = np.empty_like(target_dict['temp'])
        for i, temp in enumerate(target_dict['temp']):
            if temp == 0.0:
                target_dict['beta'][i] = 1.0/300.0
            else:
                target_dict['beta'][i] = 1.0/temp

        print(di, 'ene', target_dict['energy'])

        # Collect energy and force statistics from reference configurations
        stats_dict = {'energy':[], 'forces':[]}
        for xyz, box in zip(target_dict['xyz'], target_dict['box']):

            # calculate pair distance matrices (absalute values, components)
            if 0.5*box > sc[-1]:
                rr, rx = pair_dist(xyz, box)
            else:
                print('box', 0.5*box, sc[-1])
                rr, rx = pair_dist_minibox(xyz, box)

            print('mindist', np.where(rr > 0.0, rr, 10000.0).min())

            # calculate sufficient statistics for energies and forces from pair distances
            a1, ar, a2, f1, fr, f2 = get_stats_EAM(rr, rx, sc)
            print(xyz.shape, box)
            #print(xyz)
            print('x', a1.shape, rr.shape, np.sum(np.abs(a1)))

            stats_dict['energy'].append(np.array([ar, a2, a1]))
            stats_dict['forces'].append(np.array([fr, f2, f1]))
            #print(stats_dict['energy'][-1])

        # add dataset
        stats_data['dset'+str(di)] = stats_dict
        target_data['dset'+str(di)] = target_dict

    # pickle stats and target data to be used for optimization
    with open('../../../data/working/target.pickle', 'wb') as fo:
        pickle.dump(target_data, fo)

    with open('../../../data/working/stats.pickle', 'wb') as fo:
        pickle.dump(stats_data, fo)

