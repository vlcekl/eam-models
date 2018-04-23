#!/usr/bin/python
#
# File name:   harmdet.py
# Date:        2012/01/26 17:49
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
#import string as s
import re
import numpy as np
from scipy.optimize import fmin

temp = 298.15
beta = 1000.0/(8.314472*temp)
pref = 1.0 # reference pressure
fac = (1.602176487e-19)**2*6.02214179e23*1e10/(4.*np.pi*8.854187817e-12)/1000.0

nmax = 190000
#nmax = 2000
fn = 1.0/float(nmax)
lmin = 0
lread = 149
lmax = 130

sa = np.zeros((nmax), dtype=float)
sc = np.zeros((nmax), dtype=float)
sq = np.zeros((nmax), dtype=float)
pa = np.zeros((nmax), dtype=float)
pc = np.zeros((nmax), dtype=float)
pi = np.zeros((nmax), dtype=float)
pq = np.zeros((nmax), dtype=float)
nx = np.zeros((nmax), dtype=float)
hhh = np.zeros((lread,nmax,3), dtype=float)
hrf = np.zeros((lread,3), dtype=float)
hrsc = np.zeros((lread,3), dtype=float)
hee = np.zeros(nmax)
uuu = np.zeros((nmax), dtype=float)
eee = np.zeros((nmax), dtype=float)
dh = np.zeros((8), dtype=float)

def f_his(p, sa, sc, sq, hhh, hrf, hee, pi, pa, pc, pq, s):
    cc = 4.0*p[1]*p[0]**6
    aa = cc*p[0]**6
    qq = p[2]**2

    uuu[:] = beta*(aa*sa[:] - cc*sc[:] + qq*sq[:] - hee[:])
    uave = np.sum(uuu[:])*fn
    uuu[:] = uuu[:] - uave
    eee = np.exp(-uuu[:])
    fx = 1.0/np.sum(eee)

    # pressure loss function
    nx = pi + aa*pa/4.184 - cc*pc/4.184 + qq*pq
    na = np.sum(nx*eee)*fx
    nv = np.sum(nx*nx*eee)*fx - na*na
    s2p = np.arccos(np.exp(-0.25*(na-npref)**2/(2.0*nv)))**2
    #print *, npref, na, nv, nv**0.5, fx
    #print *, nx(ncom), pi(ncom), pa(ncom), aa/4.184, pc(ncom), cc/4.184, pq(ncom), qq

    for l in range(lmin, lmax):
        hrsc[l,0] = np.sum(hhh[l,:,0]*eee[:])*fx
        hrsc[l,1] = np.sum(hhh[l,:,1]*eee[:])*fx
        hrsc[l,2] = np.sum(hhh[l,:,2]*eee[:])*fx
    dh[0] = np.arccos(np.sum(np.sqrt(hrsc[0:lmax,0])*np.sqrt(hrf[0:lmax,0])))
    dh[1] = np.arccos(np.sum(np.sqrt(hrsc[0:lmax,1])*np.sqrt(hrf[0:lmax,1])))
    dh[2] = np.arccos(np.sum(np.sqrt(hrsc[0:lmax,2])*np.sqrt(hrf[0:lmax,2])))
    s2g = (dh[0]**2 + 1.0*dh[1]**2 + 1.0*dh[2]**2)/3.0

    fx = -np.log(fn/fx)
    eee = np.exp(0.5*(fx - uuu))
    db = -2.0*np.log(np.sum(eee)*fn)
    ge = (fx + uave)/beta

    # ab initio loss function
    fnx = 1.0/float(len(s[:,1]))
    uvdw = aa*s[:,1] - cc*s[:,2]
    ucoul = fac*(qq*s[:,10] - 0.5*qq*s[:,11] + 0.25*qq*s[:,5])
    uuh = beta*(uvdw + ucoul - s[:,0])
    uave = np.sum(uuh)*fnx
    uuh = uuh - uave
    eee = np.exp(-uuh)
    geh = -np.log(np.sum(eee)*fnx)
    eee = np.exp(0.5*(geh - uuh))
    s2a = np.arccos(np.sum(eee)*fnx)
    #db = -2.0*np.log(np.sum(eee)*fnx)
    #s2p = 0.0
    #s2g = 0.0

    s2 = (50.0*s2g + 50.0*s2p + s2a)/101.0

    print 'b', s2, p[0],p[1],p[2], db, ge, s2a, s2p, s2g, dh[0],dh[1],dh[2]
    return s2

# read parameters
fi = open(sys.argv[4], 'r')
pars = []
while 1:
    line = fi.readline()
    if not line: break
    pars.append(np.array(map(float, re.findall('\S+', line))))

pars = np.array(pars)
print 'pars.shape', pars.shape
fi.close()

# read data
fi = open(sys.argv[1], 'r')

(siref, epref, qref, npref) = map(float, re.findall('\S+', fi.readline())[1:])

for it in range(nmax):
    line = fi.readline()
    #(hee[it], sa[it], sc[it], sq[it]) = map(float, re.findall('\S+', fi.readline()))
    farr = map(float, re.findall('\S+', fi.readline()))
    hee[it] = farr[1] + farr[2]
    sa[it] = farr[3]
    sc[it] = farr[4]
    sq[it] = farr[6]
    pa[it] = farr[8]
    pc[it] = farr[9]
    pi[it] = farr[10]
    pq[it] = farr[11]
    for l in range(lread):
        hhh[l,it,:] = np.array(map(float, re.findall('\S+', fi.readline())))
    hhh[:lmax,it,0] = hhh[:lmax,it,0]/np.sum(hhh[0:lmax,it,0])
    hhh[:lmax,it,1] = hhh[:lmax,it,1]/np.sum(hhh[0:lmax,it,1])
    hhh[:lmax,it,2] = hhh[:lmax,it,2]/np.sum(hhh[0:lmax,it,2])
    print 'data', it

hee = hee*4.184
sq = sq*4.184
print 'shapes', sa.shape, sc.shape, sq.shape, hhh.shape
fi.close()

# read target href
fi = open(sys.argv[2], 'r')
for l in range(lread):
    hrf[l,:] = np.array(map(float, re.findall('\S+', fi.readline())[1:]))
    print 'href', l
hrf[:lmax,0] = hrf[:lmax,0]/sum(hrf[:lmax,0])
hrf[:lmax,1] = hrf[:lmax,1]/sum(hrf[:lmax,1])
hrf[:lmax,2] = hrf[:lmax,2]/sum(hrf[:lmax,2])
print 'shapes', hrf.shape
fi.close()

# read target ab initio
fi = open(sys.argv[3], 'r')
scl = []
while 1:
    line = fi.readline()
    if not line: break
    scl.append(np.array(map(float, re.findall('\S+', line)[1:])))
scl = np.array(scl)
print 'scl.shape', scl.shape
fi.close()

for i in range(pars.shape[0]):
    so = pars[i,0]
    eo = pars[i,1]*4.184
    qo = pars[i,2]
    par_in = [so, eo, qo]

    output = fmin(f_his, par_in, args=(sa,sc,sq,hhh,hrf,hee,pi,pa,pc,pq,scl), maxiter=10000, maxfun=10000, disp=0, full_output=1)
    xopt = output[0] 
    ofunc = output[1]
    print '# fmin ', xopt[:], ofunc, f_his(xopt[:4],sa,sc,sq,hhh,hrf,hee,pi,pa,pc,pq,scl)

