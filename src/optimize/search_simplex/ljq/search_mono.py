#!/usr/bin/python

import sys
import string as s
import re
import numpy as np
from scipy.optimize import fmin, anneal

#temp = 298.15
temp = 298.
beta = 1000.0/(8.314472*temp)
fac = (1.602176487e-19)**2*6.02214179e23*1e10/(4.*np.pi*8.854187817e-12)/1000.0

def f_db(p, s):
    """Free energy perturbation"""
    fn = 1.0/float(len(s[:,1]))

    cc = 4.0*p[1]*p[0]**6
    aa = cc*p[0]**6
    qh = p[2]
    qo = -2.0*qh

    uvdw = aa*s[:,1] - cc*s[:,2]
    ucoul = fac*(qo*qo*s[:,10] + qo*qh*s[:,11] + qh*qh*s[:,5])
    uuu = beta*(uvdw + ucoul - s[:,0])
    uave = np.sum(uuu)*fn
    uuu = uuu - uave
    eee = np.exp(-uuu)
    ge = -np.log(np.sum(eee)*fn)

    eee = np.exp(0.5*(ge - uuu))
    db = -2.0*np.log(np.sum(eee)*fn)

    print p[0], p[1], p[2], db, uave, np.sum(np.exp(-uuu))*fn, uuu[1]+uave
    return db

if __name__ == "__main__":

    fi = open(sys.argv[1], 'r')
    scl = []
    while 1:
        line = fi.readline()
        if not line: break
        scl.append(np.array(map(float, re.findall('\S+', line)[1:])))
    fi.close()
    scl = np.array(scl)
    print 'scl.shape', scl.shape

    # param values
    fi = open(sys.argv[2], 'r')
    data = []
    while 1:
        line = fi.readline()
        if not line: break
        data.append(np.array(map(float, re.findall('\S+', line))))
    fi.close()
    data = np.array(data)
    print 'data.shape', data.shape

    for i in range(data.shape[0]):
        par_in = [data[i,0], data[i,1], data[i,2]]
        output = fmin(f_db, par_in, args=(scl,), maxiter=10000, maxfun=10000, disp=0, full_output=1)
        print '# fmin ', output[0][0], output[0][1], output[0][2], output[1], f_db(xopt[:3], scl)
