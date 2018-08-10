#!/usr/bin/python

import sys
import string as s
import re
import numpy as np

if __name__ == "__main__":

    doh = 1.5
    dot = 2.7
    f = open(sys.argv[1], 'r')

    line = f.readline()
    line = f.readline()
    (ax, ay, az) = map(float, re.findall('\S+', f.readline()))
    (bx, by, bz) = map(float, re.findall('\S+', f.readline()))
    (cx, cy, cz) = map(float, re.findall('\S+', f.readline()))

    aa = (ax*ax + ay*ay + az*az)**0.5
    bb = (bx*bx + by*by + bz*bz)**0.5
    cc = (cx*cx + cy*cy + cz*cz)**0.5

    cosa = (bx*cx + by*cy + bz*cz)/(bb*cc)  # cos(alpha) angle between b and c
    cosb = (cx*ax + cy*ay + cz*az)/(cc*aa)  # cos(beta)  angle between a and c
    cosc = (ax*bx + ay*by + az*bz)/(aa*bb)  # cos(gamma) angle between a and b
    sina = (1.0 - cosa*cosa)**0.5
    sinb = (1.0 - cosb*cosb)**0.5
    sinc = (1.0 - cosc*cosc)**0.5

    lx = aa
    xy = bb*cosc
    xz = cc*cosb
    ly = (bb*bb - xy*xy)**0.5
    yz = (bb*cc*cosa - xy*xz)/ly
    lz = (cc*cc - xz*xz - yz*yz)**0.5


    # components
    axbx = (ay*bz - az*by)
    axby = (az*bx - ax*bz)
    axbz = (ax*by - ay*bx)
    nrm =  (axbx**2 + axby**2 + axbz**2)**0.5
    axbx = axbx/nrm
    axby = axby/nrm
    axbz = axbz/nrm

    dx = (axby*az - axbz*ay)
    dy = (axbz*ax - axbx*az)
    dz = (axbx*ay - axby*ax)
    nrm = (dx*dx + dy*dy + dz*dz)**0.5
    dx = dx/nrm
    dy = dy/nrm
    dz = dz/nrm

    ccx = (cx*ax + cy*ay + cz*az)/aa
    ccy = (cx*dx + cy*dy + cz*dz)
    ccz = (cx*axbx + cy*axby + cz*axbz)

    bbx = bb*cosc
    bby = bb*sinc

    aax = aa

#    print (axbx**2 + axby**2 + axbz**2)**0.5, axbx, axby, axbz
#    print (dx*dx + dy*dy + dz*dz)**0.5, dx, dy, dz
#    print 'cc', cc, (ccx*ccx + ccy*ccy + ccz*ccz)**0.5
#    print 'cx', ccx, cc*cosb

    line = f.readline()
    (npc, nph, npo, npt, npx, npw) = map(int, re.findall('\S+', f.readline()))
    line = f.readline()
    ntot = npc + nph + npo + npt + npx + npw
    #print ntot

    nbnd = 2*npw + nph + 3*npo
    nang = npw + 3*nph


    ai = np.zeros((ntot), dtype=float)
    bi = np.zeros((ntot), dtype=float)
    ci = np.zeros((ntot), dtype=float)


    nam = []
    mol = []
    for ip in range(ntot):
        (ai[ip], bi[ip], ci[ip]) = map(float, re.findall('\S+', f.readline()))
        if ip < npc:
            nam.append('C')
            mol.append(0)
        elif ip < npc+nph:
            nam.append('H')
            mol.append(0)
        elif ip < npc+nph+npo:
            nam.append('O')
            mol.append(0)
        elif ip < npc+nph+npo+npt:
            nam.append('Ti')
            mol.append(0)
        elif ip < npc+nph+npo+npt+npx:
            nam.append('Hw')
            mol.append(0)
        else:
            nam.append('Ow')
            mol.append(0)

    # find bonds and bring constrain-bonded atoms together
    ib = 0
    ia = 0
    bnd = []
    ang = []
    for ip in range(npc+nph,npc+nph+npo):  # hydroxyl oxygens
        nbh = 0
        nbt = 0
        tb = []
        hb = []
        for jp in range(npc, npc+nph+npo+npt):
            da = ai[jp] - ai[ip]
            db = bi[jp] - bi[ip]
            dc = ci[jp] - ci[ip]

            if da < -0.5:
                da = da + 1.0
            elif da > 0.5:
                da = da - 1.0

            if db < -0.5:
                db = db + 1.0
            elif db > 0.5:
                db = db - 1.0

            if dc < -0.5:
                dc = dc + 1.0
            elif dc > 0.5:
                dc = dc - 1.0

            x = da*ax + db*bx + dc*cx
            y = da*ay + db*by + dc*cy
            z = da*az + db*bz + dc*cz

            rr = x*x + y*y + z*z

            if nam[jp] == 'H' and rr < doh*doh:
               ib = ib + 1
               bnd.append([ib, 1, ip+1, jp+1])
               ai[jp] = ai[ip] + da
               bi[jp] = bi[ip] + db
               ci[jp] = ci[ip] + dc
               nbh = nbh + 1
               hb.append(jp)
#               print 'oh', rr**0.5
            elif nam[jp] == 'Ti' or nam[jp] == 'Tx':
               if rr < dot*dot:
                   ib = ib + 1
                   bnd.append([ib, 2, ip+1, jp+1])
                   ai[jp] = ai[ip] + da
                   bi[jp] = bi[ip] + db
                   ci[jp] = ci[ip] + dc
                   nam[jp] = 'Tx'
                   nbt = nbt + 1
                   tb.append(jp)
#                   print 'ot', rr**0.5

#        print ip, nbh, nbt
        ia = ia + 1
        ang.append([ia, 1, hb[0]+1, ip+1, tb[0]+1])
        ia = ia + 1
        ang.append([ia, 1, hb[0]+1, ip+1, tb[1]+1])
        ia = ia + 1
        ang.append([ia, 1, hb[0]+1, ip+1, tb[2]+1])

    im = 1
    for ip in range(npc+nph+npo+npt+npx,ntot):  # water oxygens
        nbh = 0
        hb = []
        im = im + 1
        mol[ip] = im
        for jp in range(npc+nph+npo+npt,npc+nph+npo+npt+npx):  # water hydrogens
            da = ai[jp] - ai[ip]
            db = bi[jp] - bi[ip]
            dc = ci[jp] - ci[ip]

            if da < -0.5:
                da = da + 1.0
            elif da > 0.5:
                da = da - 1.0

            if db < -0.5:
                db = db + 1.0
            elif db > 0.5:
                db = db - 1.0

            if dc < -0.5:
                dc = dc + 1.0
            elif dc > 0.5:
                dc = dc - 1.0

            x = da*ax + db*bx + dc*cx
            y = da*ay + db*by + dc*cy
            z = da*az + db*bz + dc*cz

            rr = x*x + y*y + z*z

            if nam[jp] == 'Hw' and rr < doh*doh:
               ib = ib + 1
               bnd.append([ib, 3, ip+1, jp+1])
               ai[jp] = ai[ip] + da
               bi[jp] = bi[ip] + db
               ci[jp] = ci[ip] + dc
               nbh = nbh + 1
               hb.append(jp)
               mol[jp] = mol[ip]
#               print 'oh', rr**0.5

#        print ip, nbh
        ia = ia + 1
        ang.append([ia, 2, hb[0]+1, ip+1, hb[1]+1])

    print '# MXene + H2O'
    print ''
    print '   ', ntot, 'atoms'
    print '   ', nbnd, 'bonds'
    print '   ', nang, 'angles'
    print ''
    print '   ', 7, 'atom types'
    print '   ', 3, 'bond types'
    print '   ', 2, 'angle types'
    print ''

    print '   ', 0.0, lx, 'xlo xhi'
    print '   ', 0.0, ly, 'ylo yhi'
    print '   ', 0.0, lz, 'zlo zhi'
    print '   ', xy, xz, yz, 'xy xz yz'
    print ''
    print ' Masses'
    print '  '
    print '           1   47.867'
    print '           2   47.867'
    print '           3   12.0107'
    print '           4   15.9994'
    print '           5   1.0079'
    print '           6   15.9994'
    print '           7   1.0079'
    print '          '
    print ' Pair Coeffs'
    print '          '
    print '           1 0.000000 3.166000'
    print '           2 0.000000 3.166000'
    print '           3 0.000000 3.166000'
    print '           4 0.155354 3.166000'
    print '           5 0.000000 3.166000'
    print '           6 0.155354 3.166000'
    print '           7 0.000000 3.166000'
    print '          '
    print ' Bond Coeffs'
    print '          '
    print '           1   400.00000000000        0.975'
    print '           2   400.0000000000        2.1155'
    print '           3   400.000000000000        1.0'
    print '          '
    print ' Angle Coeffs'
    print '          '
    print '           1   10.000000000000        109.47'
    print '           2   50.000000000000        109.47'


    print ''
    print 'Atoms'
    print ''

    im = 1
    for ip in range(ntot):
        x = aax*ai[ip] + bbx*bi[ip] + ccx*ci[ip]
        y =              bby*bi[ip] + ccy*ci[ip]
        z =                           ccz*ci[ip]
        #print 'xx', aax*ai[ip], bbx*bi[ip], ccx*ci[ip]
        
        if nam[ip] == 'C':
            it = 3
            q = -0.39
            print ip+1, im, it, q, x, y, z
#            print 'C', x, y, z
        elif nam[ip] == 'H':
            it = 5
            q = 0.48
            print ip+1, im, it, q, x, y, z
#            print 'H', x, y, z
        elif nam[ip] == 'O':
            it = 4
            q = -0.79
            print ip+1, im, it, q, x, y, z
#            print 'O', x, y, z
        elif nam[ip] == 'Ti':
            it = 2
            q = 0.56
            print ip+1, im, it, q, x, y, z
#            print 'Ti', x, y, z
        elif nam[ip] == 'Tx':
            it = 1
            q = 0.42
            print ip+1, im, it, q, x, y, z
#            print 'Tx', x, y, z
        elif nam[ip] == 'Ow':
            it = 6
            im = mol[ip]
            q = -0.8476
            print ip+1, im, it, q, x, y, z
#            print 'Ow', x, y, z
        elif nam[ip] == 'Hw':
            it = 7
            im = mol[ip]
            q = 0.4238
            print ip+1, im, it, q, x, y, z
#            print 'Hw', x, y, z

    print ''
    print 'Bonds'
    print ''
    for bo in bnd:
        print bo[0], bo[1], bo[2], bo[3]

    print ''
    print 'Angles'
    print ''
    for bo in ang:
        print bo[0], bo[1], bo[2], bo[3], bo[4]

    f.close()

