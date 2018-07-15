#!//anaconda/envs/py36/bin/python
#
# File name:   mkeam.py
# Date:        2018/06/27 16:42
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np

def read_params(file_name):
    """ read parameter file """

    with open(file_name, 'r') as f:
        for line in iter(f.readline, ''):
            sarr = re.findall('\S+', line)

def fembed(p, nd, dd):
    return tab

def edens(p, nr, dr):
    return tab

def upair(p, nr, dr):
    return tab

def write_lammps(fembed_tab, edens_tab, upair_tab):
    return

write_yuri(fembed_tab, edens_tab, upair_tab):
    return

if __name__ == "__main__":

    params = read_params(sys.argv[1])


    # generate embedding function table
    fembed_tab = fembed(p_fembed, n_rho, d_rho)

    # generate electronic density table
    edens_tab = edens(p_edens, n_r, d_r)

    # generate pair potential table
    upair_tab = upair(p_upair, n_r, d_r)

    # write tabulated potential
    with open(sys.argv[2], 'w') as f:
        # write
        for line in iter(f.readline, ''):
            sarr = re.findall('\S+', line)

# end of mkeam.py 
