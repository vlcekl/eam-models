#!//anaconda/envs/py36/bin/python
#
# File name:   test_pickle.py
# Date:        2019/01/10 14:49
# Author:      Lukas Vlcek
#
# Description: 
#



import sys
import re
import numpy as np
import pickle
import pprint

if __name__ == "__main__":

    with open('marinica_params.pickle', 'rb') as f:
        params = pickle.load(f)

    pprint.pprint(params)
