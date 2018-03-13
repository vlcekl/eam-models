#!/usr/bin/python
#
# File name:   sdopt.py
# Date:        2018/03/12 23:45
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

    for line in f:
        sarr = re.findall('\S+', line)

    f.close()

# end of sdopt.py 
