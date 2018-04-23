#!/usr/local/Python-2.7.6/bin/python

#from __future__ import print_function
import sys
import re
import numpy as np
from scipy import optimize

def sd_hist(p, q, grs, grp, hrs, hrp, hru):
    """Statistical distance between histograms (surface, profile)"""

    # apply bounds on parametes
    p[3] = 0.0
    #p[1] = 0.0
    p = np.where(p < -1.0, -1.0, p)
    p = np.where(p >  1.0,  1.0, p)

    # energy diference: bulk(1,2), surface(1,2), surface(1,1)
    uuu = beta*np.sum(hru*(p - q), axis=1)
    uave = np.sum(uuu)*fn
    #print(uave)
    uuu -= uave
    eee = np.exp(-uuu)
    fx = 1/np.sum(eee)
    #print('uu', np.sum(uuu), uave, np.sum(eee), fx)

    # statistical distance for surface configuration histogram
    dloss = np.arccos(np.sum(np.sqrt(np.sum(hrs*eee, axis=1)*fx*grs)))**2

    # statistical distance for concentration profile
    xp = np.sum(hrp*eee, axis=1)*fx # average concentration across layers
    vp = np.var(hrp)                # approximate unknown target variance by that of the reference 
    dloss += np.sum(np.arccos(np.exp(-(xp - grp)**2/(4.0*vp)))**2)

    fx = -np.log(fn/fx)
    eee = np.exp(0.5*(fx - uuu))
    db = -2.0*np.log(np.sum(eee)*fn)
    ge = (fx + uave)/beta

    #print('b', dloss, p[0],p[1],p[2], db, ge)
    return dloss

def dg(p, q, hru):
    """Free energy and Bhattacharyya distance between optimal and reference systems"""
    # apply bounds on parametes
    #p[2] = 0.0
    p = np.where(p < -1.0, -1.0, p)
    p = np.where(p >  1.0,  1.0, p)
    uuu = beta*np.sum(hru*(p - q), axis=1)
    uave = np.sum(uuu)*fn
    uuu -= uave
    eee = np.exp(-uuu)
    fx = -np.log(np.sum(eee)*fn)
    eee = np.exp(0.5*(fx - uuu))
    db = -2.0*np.log(np.sum(eee)*fn)
    ge = (fx + uave)/beta
    return (db, ge)

if __name__ == "__main__":

    # read target data (histograms)
    with open(sys.argv[1], 'r') as fi:
        # profile histogram
        fi.readline()
        lpmax = int(re.findall('\S+', fi.readline())[1])
        grp = np.array([float(re.findall('\S+', fi.readline())[1]) for _ in range(lpmax)])
        grp = grp[1:12]
        # surface histogram
        fi.readline()
        lsmax = int(re.findall('\S+', fi.readline())[1])
        grs = np.array([float(re.findall('\S+', fi.readline())[1]) for _ in range(lsmax)])
        grs = grs/np.sum(grs) # normalize

    # read reference data (energies, histograms)
    with open(sys.argv[2], 'r') as fi:
        hrp = []
        hrs = []
        hrv = []
        hrx = []

        line = fi.readline()
        # cycle over reference histograms
        nmax = 0
        for line in iter(fi.readline, "ENDHST\n"):
            nmax = nmax + 1
            # profile histograms to compare with concentration profile
            line = fi.readline()
            lpmax = 18 #20 #18
            hrp.append([float(re.findall('\S+', fi.readline())[2]) for _ in range(lpmax)])

            # surface histogram (now unused)
            line = fi.readline()
            lvmax = 32 
            hrv.append([float(re.findall('\S+', fi.readline())[1]) for _ in range(lvmax)])

            # surface histogram
            line = fi.readline()
            hrs.append([float(re.findall('\S+', fi.readline())[1]) for _ in range(lsmax)])

            # pair interaction histogram - bulk[0] and surface[1], 4,5, and 7 are useful
            line = fi.readline()
            lumax = 10
            hrx.append([list(map(float, re.findall('\S+', fi.readline())[2:4])) for _ in range(lumax)])

        hrp = np.array(hrp).transpose()/64.0 #64.0 256.0
        hrp = hrp[1:grp.shape[0]+1,:] # cutoff deep layers
        hrv = np.array(hrv).transpose()
        hrv = hrv*float(nmax)/np.sum(hrv)
        hrs = np.array(hrs).transpose()
        hrs = hrs*float(nmax)/np.sum(hrs)
        hrx = np.array(hrx)

    # read reference and starting parameters
    with open(sys.argv[3], 'r') as fi:
        # reference parameters
        pref = np.array(list(map(float, re.findall('\S+', fi.readline()))))
        print('pref', pref)
        # initial parametes for search
        pars = []
        for line in iter(fi.readline, ''):
            pars.append(np.array(list(map(float, re.findall('\S+', line)))))

    # chose appropriate histograms for selected parameters

    print('shapes', hrx.shape, pref.shape, len(pars))
    hru = np.zeros((hrx.shape[0],pref.shape[0]), dtype=float)
    hru[:,0] = hrx[:,5,0]
    hru[:,1] = hrx[:,5,1]
    hru[:,2] = hrx[:,4,1]
    hru[:,3] = hrx[:,8,0]

    temp = 0.35
    temp = float(sys.argv[4])
    beta = 1.0/temp
    uuu = np.zeros((nmax), dtype=float)
    eee = np.zeros((nmax), dtype=float)
    fn = 1.0/float(nmax)

    for par_in in pars: # cycle over differet starting parameter sets
        print('# start ', sd_hist(par_in, pref, grs, grp, hrs, hrp, hru), par_in, dg(par_in, pref, hru))
        output = optimize.fmin(sd_hist, par_in, args=(pref, grs, grp, hrs, hrp, hru), maxiter=100000, maxfun=10000, disp=0, full_output=1)
        print('# end   ', output[1], output[0][:], dg(output[0][:4], pref, hru))
        print('')
