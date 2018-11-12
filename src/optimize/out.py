#!//anaconda/envs/py36/bin/python
#
# File name:   optscript.py
# Date:        2018/11/08 13:16
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import os
import re
import pickle
import numpy as np
from scipy.optimize import fmin

sys.path.append('../../../statmechlib')
from statmechlib.read_write import read_vasp
from statmechlib.preprocessing import Trajectory
from statmechlib.forcefields import sd2_loss, utot_EAM

if __name__ == "__main__":


    working = '../../data/working'

    # load target data
    with open(os.path.join(working, 'output_pars.pickle'), 'rb') as fi:
        opt_pars = pickle.load(fi)
        for o in opt_pars:
            print(o[0][0],'\n', *o[1])
#        print(opt_pars[0][0], opt_pars[1])


#    with open(os.path.join(working, 'output_pars.pickle'), 'wb') as fo:
#        pickle.dump(optimal_parameters, fo)

