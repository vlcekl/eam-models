#!/usr/bin/python
#
# File name:   harmdet.py
# Date:        2012/01/26 17:49
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np
from scipy.optimize import fmin, anneal

lread = 149
nmax = 190000
uu = np.zeros(nmax)

# read data
fi = open(sys.argv[1], 'r')
for it in range(nmax):
    line = fi.readline()
    (uu[it], faux, faux, faux) = map(float, re.findall('\S+', fi.readline()))
    for l in range(lread):
        line = fi.readline()
    print 'data', it

uu = uu*4.184
nu = uu.shape[0]

uave = np.sum(uu)/float(nu)
uvar = np.sum(uu*uu)/float(nu) - uave*uave
umin = uu.min()
umax = uu.max()

print '#', uave, uvar, uvar**0.5, umin, umax, (umin+umax)*0.5

du = 10.0
lumax = int((umax - umin)/du) + 1

uhist = np.zeros((lumax), dtype=float)

for i in range(nu):
    lu = int((uu[i] - umin)/du)
    uhist[lu] = uhist[lu] + 1.0
    #print i, uu[i]

uhist = uhist/np.sum(uhist)

for lu in range(lumax):
    print float(lu)*du + umin, uhist[lu]
