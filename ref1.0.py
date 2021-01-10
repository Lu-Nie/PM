#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 21:01:26 2020

@author: njl
"""

import numpy as np
from time import time
import healpy as hp
from multiprocessing import Pool
import tables
import sqlite3
import sys
import os
#from lsd import DB
#from lsd import bounds as lsdbounds
# various functions are defined here
import pmfuns

#import calcMedianAndResiduals as ca
import pandas as pd
import matplotlib.pyplot as plt
import spherical_to_tangential as s2t
import random
import esutil
from esutil.numpy_util import match
from astropy.io import fits
from astropy.table import Table
import matplotlib
#import calcMedianAndResiduals

input_file = './data/hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f475w-f814w_v1_st.fits'
hdu_list = fits.open(input_file, memmap=True)
evt_data = Table(hdu_list[1].data)
CRA = np.median(evt_data['RA'])
CDE = np.median(evt_data['DEC'])




table = pd.read_csv('./data/B9ComB11_1.0.csv')
rootdir = "./data/" 
refDB = rootdir + "referenceGal.csv"
ref_cri = (table['objid'] != 0)

objIDarray, medianXiArray, medianEtaArray, noOfObsArray, residualXiArray, residualEtaArray, medianResidualXiArray, medianResidualEtaArray = pmfuns.calcMedianAndResiduals2(table, CRA, CDE, ref_cri, gnomonic=True, debug=False)

