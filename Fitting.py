#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 09:58:09 2020

@author: njl
"""

import numpy as np
import tables
import matplotlib.pyplot as plt
from esutil.coords import sphdist
#from astropy.time import Time
from time import time
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plotfun as pf
import plotfun2 as pf2
import sqlite3
import sys
import os
#from lsd import DB
#from lsd import bounds as lsdbounds
import healpy as hp
from multiprocessing import Pool
import spherical_to_tangential as s2t
import random
import astropy.io.fits as fits
import glob as gb
from esutil.numpy_util import match
import esutil
from astropy.time import Time
import pmfuns
from astropy.table import Table
import pandas as pd
htm_mesh = esutil.htm.HTM(10)

root = "./data/" 
input_file = './data_per/hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f475w-f814w_v1_gst.fits'
hdu_list = fits.open(input_file, memmap=True)
evt_data = Table(hdu_list[1].data)
CRA = np.median(evt_data['RA'])
CDE = np.median(evt_data['DEC'])

#plotFlag = False
plotFlag = True
chunkNo = 0
h5fileName = "correctedStarMock" + "%d" %chunkNo
h5file = tables.open_file(root+"%s.h5" %h5fileName)
table = h5file.root.correctedStar


if(1):
    try:
        h5fileNameSDSS = "correctedStarSDSSMock" + "%d" %chunkNo
        h5fileSDSS = tables.open_file(root+"%s.h5" %h5fileNameSDSS)
        tableSDSS = h5fileSDSS.root.correctedStar
    except IOError or tables.exceptions.HDF5ExtError or tables.exceptions.NoSuchNodeError:
        tableSDSS = 0


    magMask = '(mr>10.0) & (mr<30.0)' #'(mr>12.0) & (mr<22.0)''

    table0 = table.read_where(magMask)
    obj_id0 = table0['obj_id']
    ra0 = np.float64(table0['ra'])
    #raErr0 = table0['raErr']
    dec0 = np.float64(table0['dec'])
    #decErr0 = table0['decErr']
    mjd0 = table0['mjd']
    mr0 = table0['mr']
    
    h5file.close()

    if(tableSDSS==0):
        obj_id1 = np.zeros(10) - 999
        ra1 = np.zeros(10)
        raErr1 = np.zeros(10)
        dec1 = np.zeros(10)
        decErr1 = np.zeros(10)
        mjd1 = np.zeros(10)
	print "No coverage in SDSS this chunk"	
    else:
        obj_id1 = tableSDSS.col('obj_id')
        ra1 = tableSDSS.col('ra')
        raErr1 = tableSDSS.col('raErr')
        dec1 = tableSDSS.col('dec')
        decErr1 = tableSDSS.col('decErr')
        mjd1 = tableSDSS.col('mjd')
        h5fileSDSS.close()

    uniqueStar = np.unique(obj_id0)
    #m1, m2 = match(uniqueStar, obj_id1, presorted=True)
    #uniqueStar = np.unique(uniqueStar[m1])


    csst_o = pd.read_csv('./data/all.csv')
    test_cri = (csst_o['obj_id'] == uniqueStar[0])
    obj_cri = (obj_id0 == uniqueStar[0])
    print sum(test_cri), sum(obj_cri)
    print csst_o['mjd'][test_cri]-min(csst_o['mjd'][test_cri]), mjd0[obj_cri]-min(mjd0[obj_cri])
    #print csst_o['ra'][test_cri], ra0[obj_cri]

    xitmp0, etatmp0, status = s2t.ds2tp(csst_o['ra'][test_cri], csst_o['dec'][test_cri], CRA, CDE)
    xitmp0_, etatmp0_, status = s2t.ds2tp(ra0[obj_cri], dec0[obj_cri], CRA, CDE)
    print (xitmp0-min(xitmp0))*3600000, (xitmp0_-min(xitmp0_))*3600000
    print "important:", (etatmp0-min(etatmp0))*3600000, (etatmp0_-min(etatmp0_))*3600000



    #print raErr0, decErr0
    obj_id2 = np.zeros(10) - 999
    ra2 = np.zeros(10)
    raErr2 = np.zeros(10)
    dec2 = np.zeros(10)
    decErr2 = np.zeros(10)
    mjd2 = np.zeros(10)
    #rand = random.sample(range(len(uniqueStar)), 30000)
    packParameterList = [(index, uniqueID, obj_id0, ra0, \
	dec0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, dec2, decErr2, mjd2, \
	plotFlag, 0, CRA, CDE, root) for index, uniqueID in enumerate(uniqueStar[0:20])]

    pool = Pool(processes=20)
    iterator = pool.imap_unordered(pf.fittingPM3,packParameterList, chunksize=200)

    counterForRows = 0
    pm = []
    for res in iterator:
        if(len(res)==0):
            print "there is no enough observation:", res
        else:
            uniqueID_, rMagtmp_, R2XI_, chi2XI_, muXI_, muErrXI_, muXxi_, muErrXxi_, R2ETA_, \
    	        chi2ETA_, muETA_, muErrETA_, muXeta_, muErrXeta_, xpm_, flag_, ra0_, dec0_ = res
            pm.append([uniqueID_, rMagtmp_, R2XI_, chi2XI_, muXI_, muErrXI_, muXxi_, muErrXxi_, R2ETA_, \
                chi2ETA_, muETA_, muErrETA_, muXeta_, muErrXeta_, xpm_, flag_, ra0_, dec0_])
            counterForRows += 1
	#print counterForRows

    #terminate the pool of multiprocessors
    pool.terminate()

    pm = np.array(pm,dtype=object)
    #if(not plotFlag): np.save(root+"figs/PM_radec_Mock0.npy", pm)
    print "counter=total no of rows added to file", counterForRows