#!//anaconda/envs/py36/bin/python
#
# File name:   search.py
# Date:        2018/07/24 16:16
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np

def utot_EAM(params, stats):
    """
    Calculates configurational energy from EAM sufficient statistics and model parameters

    Parameters
    ----------
    params : list of lists and floats
             EAM interaction parameters (spline coefficients array and embedding function parameters)
    stats  : list of lists and floats
             Sufficient statistics for a trajectory of configurations

    Returns
    -------
    u_total: float
             total configurational energy (sum of pair and manybody interactions) for trajectory of configurations
    """

    n_sample = stats.shape[0]

    # pair interactions from array of spline coefficeints and corresponding statistic
    u_pair = np.array([sum([a*s for a, s in zip(params[0], stats[i, 0])]) for i in range(n_sample)])

    # manybody interactions from embedding function parameters and corresponding statistics
    u_many = np.array([params[1]*stats[i, 1] + params[2]*stats[i,2] for i in range(n_sample)])

    u_total = u_pair + u_many

    return u_total

def ftot_EAM(params, stats):
    """
    Calculates configurational energy from EAM sufficient statistics and model parameters

    Parameters
    ----------
    params : list of lists and floats
             EAM interaction parameters (spline coefficients array and embedding function parameters)
    stats  : list of lists and floats
             Sufficient statistics

    Returns
    -------
    f_total: float
             total configurational energy (sum of pair and manybody interactions)
    """

    # number of samples and atoms
    n_sample = stats.shape[0]
    n_atom = stats.shape[1]

    f_total = np.zeros((n_sample, 6*n_atom + 1), dtype=float)

    # cycle over samples
    for i in range(n_sample):

        # pair interactions from array of spline coefficeints and corresponding statistic
        f_pair = sum([p*s for p, s in zip(params[0], stats[i,0])]) 

        # manybody interactions from embedding function parameters and corresponding statistics
        f_many = params[1]*stats[i,1] + params[2]*stats[i,2]

        # Create a 6N + 1 array of 0, f, and -f
        f_total[1:3*n_atom+1] = f_pair.flatten() + f_many.flatten()
        f_total[3*n_atom+1:] = -f_total[1:3*n_atom+1]
            
    return f_total


def sd2_loss(params, stats, targets):
    """
    Calculates squared statistical distance loss function for configurational energies and forces.

    Parameters
    ----------
    params : list of lists and floats
             EAM interaction parameters (spline coefficients array and embedding function parameters)
    stats  : list of lists and floats
             Sufficient statistics
    targets: list of lists and floats
             target energies and forces

    Returns
    -------
    sd2: float
        squared statistical distance between model and target
    """

    # apply bounds on parametes
    #p = np.where(p < -1.0, -1.0, p)
    #p = np.where(p >  1.0,  1.0, p)

    sd2 = 0.0
    sd2f = 0.0

    # cycle over target system trajectories and statistics
    for targ, stat in zip(targets, stats):

        beta = targ['beta'] # system inverse temperature
        u_targ = targ['energy'] # target energies
        u_stat = stat['energy'] # energy statistics
        u_pars = params[0]
        n_sample = u_targ.shape[0]

        # energy diference array for a given target trajectory
        uuu = beta*(utot_EAM(u_pars, u_stat) - u_targ) # array(n_sample)
        uuu -= np.mean(uuu)
        eee = np.exp(-uuu)
 
        # are forces available?
        if 'forces' not in targ:
            
            # energy-based free energy difference and statistical distance
            ge = -np.log(np.mean(eee))   # free energy difference (shifted)
            cb = np.mean(np.exp(-0.5*(uuu - ge))) # Bhattacharyya coefficient
            sd2 += np.arccos(cb)**2              # statistical distance

        else:
            
            betad = beta*0.01  # beta * dl
            f_targ = targ['forces'] # target forces (n_sample, 1+6N) (0, 3Nf, -3Nf)
            f_stat = stat['forces'] # force statistics (n_sample, npars, 3N)
            f_pars = params[1]
        
            eeh = np.exp(-0.5*uuu)
            fff = ftot_EAM(f_pars, f_stat) # n_sample *(6N + 1) force contributions
            
            # target and model force terms
            fpave = np.mean(np.exp(betad*f_targ))
            fqave = np.mean([eee[i]*np.mean(np.exp(betad*fff[i])) for i in range(n_sample)])
            fhave = np.mean([eeh[i]*np.mean(np.exp(0.5*betad*(fff[i]+f_targ[i]))) for i in range(n_sample)])
            
            # force-based free energy difference and statistical distance
            gef = -np.log(fqave/fpave)
            cb = fhave/(fqave*fpave)**0.5
            if cb > 1: cb = 1
            sd2f += np.arccos(cb)**2

    return sd2 + sd2f


if __name__ == "__main__":

    ene_stats = np.load('ene_stats.npy')
    force_stats = np.load('force_stats.npy')

    print(ene_stats)

    print(np.sum(np.abs(force_stats)))

    # load params
    fi = open(sys.argv[2], 'r')
    para = []
    while 1:
        line = fi.readline()
        if not line: break
        para.append(np.array(map(float, re.findall('\S+', line))))

    fi.close()
    para = np.array(para)
    print 'para.shape', para.shape, para.shape[0]

    for j in range(para.shape[0]):
        # continue optimization from the grid minima
        sig = para[j,0]
        eps = para[j,1]*4.184
        qqq = para[j,2]
        kkk = para[j,3]*4.184
        rr0 = para[j,4]
        kka = para[j,5]*4.184
        rr0a = para[j,6]*3.1415926536/180.0
        par_in = [sig, eps, qqq, kkk, rr0, kka, rr0a]
        #par_in = para[j,:]
        par_0 = par_in[:]

        output = fmin(f_dhs, par_in, args=(dall_u, dall_f, par_0), maxiter=100000, maxfun=100000, disp=0, full_output=1,ftol=1e-6)
        xopt = output[0] 
        soi = xopt[0]
        eoi = xopt[1]
        ofunc = output[1]
        print '# fmin ', xopt[:7]
        print '# all ', f_dhs(xopt[:7], dall_u, dall_f, par_0)
        print '# ent ', f_ss(xopt[:7], dall_u, dall_f)

