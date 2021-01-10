#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from numpy.random import randn
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp
from esutil.numpy_util import match
from multiprocessing import Pool
import pandas as pd
import tables
#from lsd import DB
#from lsd import bounds as lsdbounds
import sys
import os  
import sqlite3
#import spherical_to_tangential as s2t
#import bootstrap as bt
from scipy.optimize import curve_fit
import esutil


# In[4]:


def calcMedianAndResiduals(table, debug):
    '''calculates median ra and dec along with residuals for a region of sky once'''
    '''table: pytables table that contains positions of PS1 objects'''
    #select rows that contain galaxies
    #maskForGal = table.get_where_list('(gal == 1)&((mjd==4000)|(mjd==3000)|(mjd==2000)|(mjd==1000))')
    maskForGal = table
    #total number of *detections* of galaxies
    noOfDetOfGal = len(table['objid'])
    noOfGal = len(table['objid'])

    #store the median ra and dec values for each galaxy 
    medianRaArray  = np.zeros(noOfGal)
    medianDecArray = np.zeros(noOfGal)

    #store the residual values for each galaxy and for each epoch
    residualRaArray  = np.zeros(noOfDetOfGal)
    residualDecArray = np.zeros(noOfDetOfGal)
    
    #store the number of observations for each galaxy 
    noOfObsArray = np.zeros(noOfGal, dtype='u1')
    
    #store the objIDs of each galaxy -- we run it in the loop just for debugging purposes 
    #needs 64-bit unsigned integer datatype
    objIDarray = np.zeros(noOfGal, dtype='i8')
    
    # pluck out objIDs, RAs, DECs, and MJDs for each detection of galaxy
    # (works *much* faster than accessing individual pytable rows)
    objIDs  = table['objid']
    #print "len(objIDs)", len(objIDs)
    #print "objIDs.size", objIDs.size
    raObjs  = table['ra']
    decObjs = table['dec']
    mjdObjs = table['mjd']

    try:
        # store the first object values in here
        currObjID = objIDs[0]
        currRaObj  = [raObjs[0]]
        currDecObj = [decObjs[0]]
        currMjdObj = [mjdObjs[0]]
        noOfObsArray[0] = 1

        #set variables
        objIDcounter = 0
        pos = 0
        t1 = time()

        #it should start from the second row, so add one.
        for i in np.arange(noOfDetOfGal-1)+1:
            
            if debug: print 'on the row number',i
            objID = objIDs[i]
            ra    = raObjs[i]
            dec   = decObjs[i]
            mjd   = mjdObjs[i]
            
            if debug: print 'object ID',objID
            if debug: print 'ra of object',ra
            if debug: print 'dec of object',dec
            if debug: print 'mjd of object',mjd
            if (objID != currObjID):
                currRaObj.append(ra)
                currDecObj.append(dec)
                currMjdObj.append(mjd)

                if debug: print 'same object now'
                noOfObsArray[objIDcounter] += 1
            else:
                if debug: print 'different object now'
                if debug: print 'run functions for the', objIDcounter, 'object'
                if(noOfObsArray[objIDcounter] >= 4.0):
                    if debug: print 'no of obs ', noOfObsArray[objIDcounter]
                    # convert list into array
                    currRaObj = np.array(currRaObj)
                    currDecObj = np.array(currDecObj)
                    
                    if debug: print 'current objID ra dec mjd :',currObjID, currRaObj, currDecObj, currMjdObj
                    # calculate and store medians
                    medianRa = np.median(currRaObj)
                    medianDec = np.median(currDecObj)

                    if debug: print 'medianRa is', medianRa
                    if debug: print 'medianDec is', medianDec
                    # store median value for each object
                    medianRaArray[objIDcounter] = medianRa
                    medianDecArray[objIDcounter] = medianDec
                    objIDarray[objIDcounter] = currObjID
                    # calculate residual
                    residualRa = currRaObj - medianRa
                    residualDec = currDecObj - medianDec
                    if debug: print 'residualRa is', residualRa
                    if debug: print 'residualDec is', residualDec
                    if debug: print 'pos is', pos
                    if debug: print 'i is ',i
                    # store residuals in separate arrays
                    residualRaArray[pos:i] = residualRa
                    residualDecArray[pos:i] = residualDec

                # store the objID no matter what
                objIDarray[objIDcounter] = currObjID
                objIDcounter += 1
                currObjID = objID
                pos = i
                # start new ra, dec, and mjd arrays and
                # don't forget to increment the noOfObsArray for this object
                currRaObj = [raObjs[i]]
                currDecObj = [decObjs[i]]
                currMjdObj = [mjdObjs[i]]
                noOfObsArray[objIDcounter] += 1

        # after going through all of the rows, make sure to process the last
        # object, if it has enough observations
        i +=1
        print 'processing the last object now'
        if(noOfObsArray[objIDcounter] >= 4.0):
            if debug: print 'obs greater than three, adding the last object too'
            currRaObj = np.array(currRaObj)
            currDecObj = np.array(currDecObj)
            # calculate and store medians
            medianRa = np.median(currRaObj)
            medianDec = np.median(currDecObj)

            if debug: print 'medianRa is', medianRa
            if debug: print 'medianDec is', medianDec
            medianRaArray[objIDcounter] = medianRa
            medianDecArray[objIDcounter] = medianDec
            # calculate residual
            residualRa = currRaObj - medianRa
            residualDec = currDecObj - medianDec
            if debug: print 'residualRa is', residualRa
            if debug: print 'residualDec is', residualDec
            if debug: print 'pos is', pos
            if debug: print 'i is ',i
            residualRaArray[pos : i] = residualRa
            residualDecArray[pos : i] = residualDec
            # also need to store the objID, no matter what
            objIDarray[objIDcounter] = currObjID
       
        print 'total time taken is', time() - t1


    except ValueError:
        print "currObjID medianRa medianDec residualRa residualDec, pos, i, objIDcounter", currObjID, medianRa, medianDec, residualRa, residualDec, pos, i, objIDcounter
        pass
        
    return objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray

# In[ ]:




