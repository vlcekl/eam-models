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
            aa[i] = sum([(rc - r)**3 for r in rr[i] if r < rc])

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

    # Lists of dictionaries for statistics and target data to be used in optimization
    stats_list = target_list = []

    # processed data directory
    target_proc = '../../../data/target_processed'

    # set of spline nodes (in parameter file?)
    #sc = [2.7, 3.252, 3.804, 4.356, 4.908, 5.46]
    sc = [2.56, 2.73, 3.252, 3.804, 4.20, 4.77]

    #files = ['marini_dft.h5', 'german_dft.h5', 'md_t300K.h5']
    files = ['marini_dft.h5']
    
    # cycle over target data and save statistics
    for file in files:

        target_dict = {}
        stats_dict = {'energy':[], 'forces':[]}

        with h5py.File(os.path.join(target_proc, file), "r") as fi:
            box = fi['box'][...]
            xyz = fi['coordinates'][...]

            target_dict['energy'] = fi['energies'][...]

            if 'forces' in fi.keys():
                # read and transform forces into (6N+1) arrays
                target_dict['forces'] = force_targ(fi['forces'][...])

                print('fshape', target_dict['forces'].shape)

        for i in range(xyz.shape[0]):

            # calculate pair distance matrices (absalute values, components)
            rr, rx = pair_dist(xyz[i], box[i])

            # calculate sufficient statistics for energies and forces from pair distances
            a1, ar, a2, f1, fr, f2 = get_stats_EAM(rr, rx, sc)
            print('x',i, a1.shape, np.sum(np.abs(a1)))

            stats_dict['energy'].append(np.array([a1, ar, a2]))
            stats_dict['forces'].append(np.array([f1, fr, f2]))

        # Add a reference configuration for zero of energy
        stats_list.append(stats_dict)
        target_list.append(target_dict)

    # pickle stats and target data to be used for optimization
    with open('stats.pickle', 'wb') as fo:
        pickle.dump(stats_list, fo)

    with open('target.pickle', 'wb') as fo:
        pickle.dump(target_list, fo)

