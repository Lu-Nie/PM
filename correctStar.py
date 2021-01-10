#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 15:54:45 2020

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
from esutil.coords import sphdist

nside = 2**24# pixelizes the sky
rootdir = "./data/" 

input_file = './data_per/hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f475w-f814w_v1_gst.fits'
hdu_list = fits.open(input_file, memmap=True)
evt_data = Table(hdu_list[1].data)
CRA = np.median(evt_data['RA'])
CDE = np.median(evt_data['DEC'])

table_gal = pd.read_csv('./data/all.csv')
table_star = pd.read_csv('./data/all.csv')

chunkNo = 0
refDB = rootdir + "referenceGalDB%d"%chunkNo
correctedStarsName = "correctedStarMock" + "%d" %chunkNo
h5correctedStarsFile = rootdir + "%s.h5" %correctedStarsName

maskForGal  = (table_gal['gal']==1)
maskForStar = (table_star['gal']==0)


galaxyIDs = table_gal['obj_id'][maskForGal]
raGalaxy   = table_gal['ra'][maskForGal]
decGalaxy  = table_gal['dec'][maskForGal]
mjdsGalaxy = table_gal['mjd'][maskForGal]


starIDs  = table_star['obj_id'][maskForStar]
raStar   = table_star['ra'][maskForStar]
decStar  = table_star['dec'][maskForStar]
#raErrStar   = table['raErr'][maskForStar]#/(3600*1000)
#decErrStar  = table['decErr'][maskForStar]#/(3600*1000)
mjdsStar = table_star['mjd'][maskForStar]
rMagStar = table_star['mr'][maskForStar]
nObsStar = np.zeros(sum(maskForStar)) + 2 #table.col('nObs')[maskForStar]


#print 'detection number for all the objs', len(table['obj_id'])
print 'detection number for all the gals', len(galaxyIDs)
print 'detection number for all the stars', len(starIDs)
print 'raStar, decStar', raStar[0:2], decStar[0:2]


mjdSortedStar  = np.sort(mjdsStar)
deltaTStar     = mjdSortedStar[0:-1] - mjdSortedStar[1:]
tBreakAtStar   = np.where(deltaTStar < (-100.0))[0]
mjdBreakAtStar = mjdSortedStar[tBreakAtStar + 1]

objIDarray, finalRaArray, finalDecArray, medianRaError, medianDecError = pmfuns.readDB(refDB)



phiForObj   = (finalRaArray*np.pi)/180 
thetaForObj = (90 - finalDecArray)* (np.pi/180) 
pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 	    
pixelIndexArray = np.unique(pixelIndexForObj)
theta, phi = hp.pix2ang(nside, pixelIndexArray)
pixelRa  = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi
print "Going to process %d pixels." % pixelIndexArray.size
phiForStar   = (raStar*np.pi)/180 
thetaForStar = (90 - decStar)* (np.pi/180) 
pixelIndexForStar = hp.ang2pix(nside, thetaForStar, phiForStar) 
print "mjdBreakAtStar:", mjdBreakAtStar
gnomonic = True

pixelIndexForStar_l = []
for i in pixelIndexForStar:
    pixelIndexForStar_l.append(int(i*10**(-9)))
print len(pixelIndexForStar_l)

packParameterList = [ (pickPixelNo, pixelRa[index], pixelDec[index],objIDarray, finalRaArray, finalDecArray, \
np.array(galaxyIDs), np.array(raGalaxy), np.array(decGalaxy),mjdsGalaxy , starIDs, mjdsStar, raStar, decStar, \
rMagStar, nObsStar, pixelIndexForStar, pixelIndexForStar_l, mjdBreakAtStar, gnomonic, CRA, CDE) \
for index, pickPixelNo in enumerate(pixelIndexArray)]





pmfuns.correctedStars(packParameterList[0:10000], h5correctedStarsFile)