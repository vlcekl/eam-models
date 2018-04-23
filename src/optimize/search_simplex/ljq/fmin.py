#!/usr/bin/python
#
# File name:   fmin.py
# Date:        2015/10/11 10:13
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np

if __name__ == "__main__":

    f = open(sys.argv[1], 'r')

    smin = []
    dmin = 10000.0
    for line in f:
        sarr = re.findall('\S+', line)
        d = float(sarr[3]) + float(sarr[4])
        #d = float(sarr[3])
        #d = float(sarr[4])
        if d < dmin:
            dmin = d
            smin = sarr
            print dmin, smin

    f.close()

    print dmin, smin

# end of fmin.py 
