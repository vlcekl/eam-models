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

uu = []
vv = []
fi = open(sys.argv[1], 'r')
it = 0
while 1:
    line = fi.readline()
    if not line: break
    uu.append(float(re.findall('\S+', line)[2]))
    vv.append(float(re.findall('\S+', line)[11]))
    it = it + 1
    print 'data', it

uu = np.array(uu)
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

vv = np.array(vv)
nv = vv.shape[0]
vave = np.sum(vv)/float(nv)
vvar = np.sum(vv*vv)/float(nv) - vave*vave
vmin = vv.min()
vmax = vv.max()
print '#', vave, vvar, vvar**0.5, vmin, vmax, (vmin+vmax)*0.5

dv = 10.0
lvmax = int((vmax - vmin)/dv) + 1

vhist = np.zeros((lvmax), dtype=float)

for i in range(nv):
    lv = int((vv[i] - vmin)/dv)
    vhist[lv] = vhist[lv] + 1.0

vhist = vhist/np.sum(vhist)
for lv in range(lvmax):
    print float(lv)*dv + vmin, vhist[lv]

