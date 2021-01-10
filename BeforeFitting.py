#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 09:30:46 2020

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

import calcMedianAndResiduals as ca
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

# the center of chunck
CRA = 11.20
CDEC = 41.62
table = pd.read_csv('./data/B9_gal.csv')
table2 = pd.read_csv('./data/B11_gal.csv')
rootdir = "./data/" 
refDB = rootdir + "referenceGalMockDB.csv"
ref_cri = (table['objid'] != 0)
ref_cri2 = (table2['objid'] != 0)

objIDarray, medianXiArray, medianEtaArray, noOfObsArray, residualXiArray, residualEtaArray, medianResidualXiArray, medianResidualEtaArray = pmfuns.calcMedianAndResiduals2(table, CRA, CDEC, ref_cri, gnomonic=True, debug=False)
objIDarray2, medianXiArray2, medianEtaArray2, noOfObsArray2, residualXiArray2, residualEtaArray2, medianResidualXiArray2, medianResidualEtaArray2 = pmfuns.calcMedianAndResiduals2(table2, CRA, CDEC, ref_cri2, gnomonic=True, debug=False)