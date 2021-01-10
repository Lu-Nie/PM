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
import spherical_to_tangential as s2t
#import bootstrap as bt
from scipy.optimize import curve_fit
import esutil
htm_mesh = esutil.htm.HTM(10)

def func(x, a, b):
    return a*x + b

def func2(x, b): # a is included in x
    return x + b
    
NGAL = 600

##################################################################
def chunksBound(raMin, decMin, raMax, decMax, chunkSize, OLdeg):
	someBigNumber = 5000
	chunkNo = 0
	#arrays to contain the bounds
	#i donot know the size of the arrays, so i specify some big number! 
	raMinArray  = np.zeros(someBigNumber)
	raMaxArray  = np.zeros(someBigNumber)
	decMinArray = np.zeros(someBigNumber)
	decMaxArray = np.zeros(someBigNumber)

	#assign the current min/max ra/dec values to the 
	currentRaMin  = raMin
	currentDecMin = decMin
	currentDecMax = decMin + chunkSize
	decAvg = (currentDecMin + currentDecMax)/2.0
	currentRaMax  = raMin + chunkSize/np.cos(np.radians(decAvg))
	raMinArray[0] = raMin
	decMinArray[0]= decMin
	raMaxArray[0] = currentRaMax
	decMaxArray[0]= currentDecMax

	#loop to calculate bounds for the chunks of sky 
	while(currentDecMax < decMax):
	    decAvg = (currentDecMin + currentDecMax)/2.0
	    print "decAvg", decAvg
	    print "cosine of decAvg", np.cos(np.radians(decAvg))
	    counter = 0 
	    while(currentRaMax < raMax):	
		chunkNo += 1
		raMinArray[chunkNo] = currentRaMax - (OLdeg/np.cos(np.radians(decAvg)))
		raMaxArray[chunkNo] = raMinArray[chunkNo] + chunkSize/np.cos(np.radians(decAvg))
		decMinArray[chunkNo]= currentDecMin
		decMaxArray[chunkNo]= currentDecMax
		counter +=1# to count after how many chunks in ra do we reach the end
		if (raMaxArray[chunkNo] > raMax):
		    print "hi", counter
		    raMaxArray[chunkNo] = raMaxArray[chunkNo] - raMax
		    raMinArray[chunkNo] = raMinArray[chunkNo] - raMax
		    break #break statement breaks out of the smallest enclosing of the while/for loop
		if (decMaxArray[chunkNo] > decMax):
		    decMaxArray[chunkNo] = decMax
		currentRaMax = raMaxArray[chunkNo]
	    currentRaMin  = raMin
	    currentRaMax  = raMin + chunkSize/np.cos(np.radians(decAvg))
	    currentDecMin = currentDecMax - OLdeg
	    currentDecMax = currentDecMin + chunkSize
	    print "chunk no", chunkNo  

	#shorten the arrays to the required length
	#add one so that it takes elements upto the chunkNo
	raMinArray  = raMinArray[0:chunkNo+1]
	raMaxArray  = raMaxArray[0:chunkNo+1]
	decMinArray = decMinArray[0:chunkNo+1]
	decMaxArray = decMaxArray[0:chunkNo+1]
	return raMinArray,raMaxArray, decMinArray, decMaxArray

##################################################################
def executeChunkWise(raMin, decMin, raMax, decMax, query, db, Star, h5file):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.

    # open a pytable file
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=80563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
	star['mr']  = row[9]# Added by HJT
        counterForRows += 1
        for j in range(len(dtype)-3): 
            star[dtype[j+1]] = row[j+1]
        if (row[7] > 0.3) & (row[7] < 1.0) & (row[8] > 0.3) & (row[8] < 1.0):
            star['gal'] = 1
        else:
            star['gal'] = 0
        star.append()
        if (counterForRows % 100000 == 0): 
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    if(noOfRowsInTable<1): 
	h5file.close()
	return noOfRowsInTable, counterForRows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

    #should return the file name
    return noOfRowsInTable, counterForRows

##################################################################
def executeChunkWiseSDSS(raMin, decMin, raMax, decMax, query, db, Star, h5file):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.

    # open a pytable file
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
        counterForRows += 1
        for j in range(len(dtype)-1): 
            star[dtype[j+1]] = row[j+1] 
        star.append()
        if (counterForRows % 100000 == 0): 
	    #print "records of SDSS stars:", row[0], row[1], row[2], row[3] , row[4], row[5], row[6], row[7]
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    if(noOfRowsInTable<1): 
	h5file.close()
	return noOfRowsInTable, counterForRows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

    #should return the file name
    return noOfRowsInTable, counterForRows

def executeChunkWise2MASS(raMin, decMin, raMax, decMax, query, db, Star, h5file):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.

    # open a pytable file
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
	star['mr']  = row[8]# Added by HJT
        counterForRows += 1
        for j in range(len(dtype)-3): 
            star[dtype[j+1]] = row[j+1]
        if (row[6] > 0.3) & (row[6] < 1.0) & (row[7] > 0.3) & (row[7] < 1.0):
            star['gal'] = 1
        else:
            star['gal'] = 0
        star.append()
        if (counterForRows % 100000 == 0): 
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    if(noOfRowsInTable<1): 
	h5file.close()
	return noOfRowsInTable, counterForRows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

    #should return the file name
    return noOfRowsInTable, counterForRows


def executeChunkWiseGAIA(raMin, decMin, raMax, decMax, query, db, Star, h5file):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.

    # open a pytable file
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
	star['mr']  = row[8]# Added by HJT
        counterForRows += 1
        for j in range(len(dtype)-3): 
            star[dtype[j+1]] = row[j+1]
        if (row[6] > 0.3) & (row[6] < 1.0) & (row[7] > 0.3) & (row[7] < 1.0):
            star['gal'] = 1
        else:
            star['gal'] = 0
        star.append()
        if (counterForRows % 100000 == 0): 
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    if(noOfRowsInTable<1): 
	h5file.close()
	return noOfRowsInTable, counterForRows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

    #should return the file name
    return noOfRowsInTable, counterForRows

def executeChunkWiseGAIA2(raMin, decMin, raMax, decMax, query, db, Star, h5file):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.

    # open a pytable file
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
	star['mr']  = row[12]# Added by HJT
        counterForRows += 1
        for j in range(len(dtype)-3): 
            star[dtype[j+1]] = row[j+1]
        if (row[10] > 0.3) & (row[10] < 1.0) & (row[11] > 0.3) & (row[11] < 1.0):
            star['gal'] = 1
        else:
            star['gal'] = 0
        star.append()
        if (counterForRows % 100000 == 0): 
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    if(noOfRowsInTable<1): 
	h5file.close()
	return noOfRowsInTable, counterForRows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

    #should return the file name
    return noOfRowsInTable, counterForRows

##################################################################
def calcMedianAndResiduals(table, debug):
    '''calculates median ra and dec along with residuals for a region of sky once'''
    '''table: pytables table that contains positions of PS1 objects'''
    #select rows that contain galaxies
    maskForGal = table.get_where_list('(gal == 1)&((mjd==4000)|(mjd==3000)|(mjd==2000)|(mjd==1000))')
    
    #total number of *detections* of galaxies
    noOfDetOfGal = maskForGal.size
    noOfGal = np.unique(table.col('obj_id')[maskForGal]).size

    #store the median ra and dec values for each galaxy 
    medianRaArray  = np.zeros(noOfGal)
    medianDecArray = np.zeros(noOfGal)

    ##################### updated by HJT#################
    medianResidualRaArray  = np.zeros(noOfGal)
    medianResidualDecArray = np.zeros(noOfGal) 
    #####################################################

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
    objIDs  = table.col('obj_id')[maskForGal]
    #print "len(objIDs)", len(objIDs)
    #print "objIDs.size", objIDs.size
    raObjs  = table.col('ra')[maskForGal]
    decObjs = table.col('dec')[maskForGal]
    mjdObjs = table.col('mjd')[maskForGal]

    ##################### updated by HJT#################
    #raErrObjs  = table.col('raErr')[maskForGal]  
    #decErrObjs = table.col('decErr')[maskForGal] 
    #####################################################

    try:
        # store the first object values in here
        currObjID = objIDs[0]
        currRaObj  = [raObjs[0]]
        currDecObj = [decObjs[0]]
	##################### updated by HJT#################
    	#currRaErrObj  = [raErrObjs[0]]
        #currDecErrObj = [decErrObjs[0]] 
    	#####################################################
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
    	    ##################### updated by HJT#################
    	    #raErr  = raErrObjs[i]  
    	    #decErr = decErrObjs[i] 
    	    #####################################################

            if debug: print 'object ID',objID
            if debug: print 'ra of object',ra
            if debug: print 'dec of object',dec
            if debug: print 'mjd of object',mjd
            if (objID == currObjID):
                currRaObj.append(ra)
                currDecObj.append(dec)
                currMjdObj.append(mjd)
    	        ##################### updated by HJT#################
    	        #currRaErrObj.append(raErr)
                #currDecErrObj.append(decErr) 
    	        #####################################################

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
		    ##################### updated by HJT#################
    	            #currRaErrObj = np.array(currRaErrObj)
                    #currDecErrObj = np.array(currDecErrObj) 
    	            #####################################################

                    if debug: print 'current objID ra dec mjd :',currObjID, currRaObj, currDecObj, currMjdObj
                    # calculate and store medians
                    medianRa = np.median(currRaObj)
                    medianDec = np.median(currDecObj)
		    ##################### updated by HJT#################
    	            #medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
                    #medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	            #####################################################

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
		    ##################### updated by HJT#################
		    # store median value for each object
                    medianResidualRaArray[objIDcounter] = np.median(residualRa)
                    medianResidualDecArray[objIDcounter] = np.median(residualDec)
		    #####################################################

                # store the objID no matter what
                objIDarray[objIDcounter] = currObjID
                objIDcounter += 1
                currObjID = objID
                pos = i
                # start new ra, dec, and mjd arrays and
                # don't forget to increment the noOfObsArray for this object
                currRaObj = [raObjs[i]]
                currDecObj = [decObjs[i]]
		##################### updated by HJT#################
    	        #currRaErrObj = [raErrObjs[i]]
                #currDecErrObj = [decErrObjs[i]]
    	        #####################################################
                currMjdObj = [mjdObjs[i]]
                noOfObsArray[objIDcounter] += 1

        # after going through all of the rows, make sure to process the last
        # object, if it has enough observations
        i +=1
        print 'processing the last object now'
        if(noOfObsArray[objIDcounter] >= 2.0):
            if debug: print 'obs greater than three, adding the last object too'
            currRaObj = np.array(currRaObj)
            currDecObj = np.array(currDecObj)
	    ##################### updated by HJT#################
    	    #currRaErrObj = np.array(currRaErrObj)
            #currDecErrObj = np.array(currDecErrObj)
    	    #####################################################
            # calculate and store medians
            medianRa = np.median(currRaObj)
            medianDec = np.median(currDecObj)

	    ##################### updated by HJT#################
    	    #medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
            #medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	    #####################################################
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
	    ##################### updated by HJT#################
	    # store median value for each object
            medianResidualRaArray[objIDcounter] = np.median(residualRa)
            medianResidualDecArray[objIDcounter] = np.median(residualDec)
	    #####################################################
            # also need to store the objID, no matter what
            objIDarray[objIDcounter] = currObjID
       
        print 'total time taken is', time() - t1


    except ValueError:
        print "currObjID medianRa medianDec residualRa residualDec, pos, i, objIDcounter", currObjID, medianRa, medianDec, residualRa, residualDec, pos, i, objIDcounter
        pass
        
    return objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray, medianResidualRaArray, medianResidualDecArray

##################################################################
def calcMedianAndResiduals2(table, ra0, dec0, ref_cri, gnomonic, debug):
    '''calculates median xi and eta along with residuals for a region of sky once'''
    '''table: pytables table that contains positions of PS1 objects'''
    #select rows that contain galaxies
    maskForGal = ref_cri#table.get_where_list(ref_cri)
    
    #total number of *detections* of galaxies
    noOfDetOfGal = sum(maskForGal)
    noOfGal = np.unique(table['obj_id'][maskForGal]).size
    #print(table['obj_id'][maskForGal])
    #print "noOfGal:", noOfGal
    #print "noOfDetOfGal:", noOfDetOfGal
    #store the median ra and dec values for each galaxy 
    medianRaArray  = np.zeros(noOfGal)
    medianDecArray = np.zeros(noOfGal)

    ##################### updated by HJT#################
    medianResidualRaArray  = np.zeros(noOfGal)
    medianResidualDecArray = np.zeros(noOfGal) 
    #####################################################

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
    objIDs  = table['obj_id'][maskForGal]
    #print "len(objIDs)", len(objIDs)
    #print "objIDs.size", objIDs.size
    raObjs  = table['ra'][maskForGal]
    decObjs = table['dec'][maskForGal]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raObjs, decObjs, status = s2t.ds2tp(raObjs, decObjs, ra0, dec0) # transform the (RA, DEC)/degree into (xi, eta)/degree, but here (xi, eta) denoted by (raObjs, decObjs)
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    mjdObjs = table['mjd'][maskForGal]

    ##################### updated by HJT#################
    #raErrObjs  = table['raErr'][maskForGal]  
    #decErrObjs = table['decErr'][maskForGal] 
    #####################################################

    try:
        # store the first object values in here
        currObjID = objIDs[0]
        currRaObj  = [raObjs[0]]
        currDecObj = [decObjs[0]]
	##################### updated by HJT#################
    	#currRaErrObj  = [raErrObjs[0]]
        #currDecErrObj = [decErrObjs[0]] 
    	#####################################################
        currMjdObj = [mjdObjs[0]]
        noOfObsArray[0] = 1

        #set variables
        objIDcounter = 0
        pos = 0
        t1 = time()

        #it should start from the second row, so add one.
        for i in np.arange(noOfDetOfGal):
            if debug: print 'on the row number',i
            objID = objIDs[i]
            ra    = raObjs[i]
            dec   = decObjs[i]
            mjd   = mjdObjs[i]
    	    ##################### updated by HJT#################
    	    #raErr  = raErrObjs[i]  
    	    #decErr = decErrObjs[i] 
    	    #####################################################

            if debug: print 'object ID',objID
            if debug: print 'ra of object',ra
            if debug: print 'dec of object',dec
            if debug: print 'mjd of object',mjd
            if (objID == currObjID):
                currRaObj.append(ra)
                currDecObj.append(dec)
                currMjdObj.append(mjd)
    	        ##################### updated by HJT#################
    	        #currRaErrObj.append(raErr)
                #currDecErrObj.append(decErr) 
    	        #####################################################

                if debug: print 'same object now'
                noOfObsArray[objIDcounter] += 1
            else:
                if debug: print 'different object now'
                if debug: print 'run functions for the', objIDcounter, 'object'
		currMjdSorted  = np.sort(currMjdObj)
    		deltaT = currMjdSorted[0:-1] - currMjdSorted[1:]
    		tBreakAt   = np.where(deltaT < (-100.0))[0]
    		currMjdBreakAt = currMjdSorted[tBreakAt+1]
		if debug: print "currMjdObj:", currMjdObj
                if(len(currMjdBreakAt) >= 1.0):#(noOfObsArray[objIDcounter] >= 2.0):
		    #print "tian test1,currMjdBreakAt", currMjdBreakAt
                    if debug: print 'no of obs ', noOfObsArray[objIDcounter]
                    # convert list into array
                    currRaObj = np.array(currRaObj)
                    currDecObj = np.array(currDecObj)
		    ##################### updated by HJT#################
		    currMjdObj = np.array(currMjdObj)
    	            #currRaErrObj = np.array(currRaErrObj)
                    #currDecErrObj = np.array(currDecErrObj)
		    '''currRaObj = currRaObj[currMjdObj<=currMjdBreakAt[2]]  ## Here notice we only choose the first 4 epochs 
		    currDecObj = currDecObj[currMjdObj<=currMjdBreakAt[2]]
		    currRaErrObj = currRaErrObj[currMjdObj<=currMjdBreakAt[2]]
		    currDecErrObj = currDecErrObj[currMjdObj<=currMjdBreakAt[2]]
		    currMjdObj = currMjdObj[currMjdObj<=currMjdBreakAt[2]]'''
    	            #####################################################

                    if debug: print 'current objID ra dec mjd :',currObjID, currRaObj, currDecObj, currMjdObj
                    # calculate and store medians
                    medianRa = np.median(currRaObj)
                    medianDec = np.median(currDecObj)
		    ##################### updated by HJT#################
    	            #medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
                    #medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	            #####################################################

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
		    indTmp = np.arange(pos,i)
		    #indTmp = indTmp[currMjdObj<=currMjdBreakAt[2]]
                    residualRaArray[indTmp] = residualRa
                    residualDecArray[indTmp] = residualDec
		    ##################### updated by HJT#################
		    # store median value for each object
                    medianResidualRaArray[objIDcounter] = np.median(residualRa)
                    medianResidualDecArray[objIDcounter] = np.median(residualDec)
		    #####################################################

                # store the objID no matter what
                objIDarray[objIDcounter] = currObjID
                objIDcounter += 1
                currObjID = objID
                pos = i
                # start new ra, dec, and mjd arrays and
                # don't forget to increment the noOfObsArray for this object
                currRaObj = [raObjs[i]]
                currDecObj = [decObjs[i]]
		##################### updated by HJT#################
    	        #currRaErrObj = [raErrObjs[i]]
                #currDecErrObj = [decErrObjs[i]]
    	        #####################################################
                currMjdObj = [mjdObjs[i]]
                noOfObsArray[objIDcounter] += 1
                #objIDcounter += 1

        # after going through all of the rows, make sure to process the last
        # object, if it has enough observations
        i +=1

        print 'processing the last object now'
        currMjdSorted  = np.sort(currMjdObj)
    	deltaT = currMjdSorted[0:-1] - currMjdSorted[1:]
    	tBreakAt   = np.where(deltaT < (-100.0))[0]
    	currMjdBreakAt = currMjdSorted[tBreakAt+1]
        if(len(currMjdBreakAt) >= 1.0):#(noOfObsArray[objIDcounter] >= 2.0):
            if debug: print 'obs greater than three, adding the last object too'
            currRaObj = np.array(currRaObj)
            currDecObj = np.array(currDecObj)
	    ##################### updated by HJT#################
    	    #currRaErrObj = np.array(currRaErrObj)
            #currDecErrObj = np.array(currDecErrObj)
	    currMjdObj = np.array(currMjdObj)
	    '''currRaObj = currRaObj[currMjdObj<=currMjdBreakAt[2]]  ## Here notice we only choose the first 4 epochs 
	    currDecObj = currDecObj[currMjdObj<=currMjdBreakAt[2]]
	    currRaErrObj = currRaErrObj[currMjdObj<=currMjdBreakAt[2]]
	    currDecErrObj = currDecErrObj[currMjdObj<=currMjdBreakAt[2]]
	    currMjdObj = currMjdObj[currMjdObj<=currMjdBreakAt[2]]'''
    	    #####################################################
            # calculate and store medians
            medianRa = np.median(currRaObj)
            medianDec = np.median(currDecObj)

	    ##################### updated by HJT#################
    	    #medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
            #medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	    #####################################################
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
            indTmp = np.arange(pos,i)
	    #indTmp = indTmp[currMjdObj<=currMjdBreakAt[2]]
            residualRaArray[indTmp] = residualRa
            residualDecArray[indTmp] = residualDec
	    ##################### updated by HJT#################
	    # store median value for each object
            medianResidualRaArray[objIDcounter] = np.median(residualRa)
            medianResidualDecArray[objIDcounter] = np.median(residualDec)
	    #####################################################
            # also need to store the objID, no matter what
            objIDarray[objIDcounter] = currObjID
       
        print 'total time taken is', time() - t1


    except ValueError:
        print "currObjID medianRa medianDec residualRa residualDec, pos, i, objIDcounter", currObjID, medianRa, medianDec, residualRa, residualDec, pos, i, objIDcounter
        pass
        
    return objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray, medianResidualRaArray, medianResidualDecArray


##################################################################
def calcMedianAndResidualsSDSS(table, tableSDSS, ra0, dec0, gnomonic, debug):
    '''calculates median xi and eta along with residuals for a region of sky once'''
    '''table: pytables table that contains positions of PS1 objects'''
    #select rows that contain galaxies
    #maskForGal = table.get_where_list('(gal == 1)&(mr<21.5)&(mr>13.5)')#('(gal == 1)&(mjd<4001)&(mjd>890)') #
    maskForGalSDSS  = tableSDSS.get_where_list('(sdss_type==3)&(mr<22.0)&(mr>13.5)')  
    m1, m2, d12 = htm_mesh.match(table.col('ra'), table.col('dec'), tableSDSS.col('ra')[maskForGalSDSS], tableSDSS.col('dec')[maskForGalSDSS], 1.5/3600, maxmatch=100000)  
    maskForGal = m1

    #total number of *detections* of galaxies
    noOfDetOfGal = maskForGal.size
    noOfGal = np.unique(table.col('obj_id')[maskForGal]).size
    #print "noOfGal:", noOfGal
    #print "noOfDetOfGal:", noOfDetOfGal
    #store the median ra and dec values for each galaxy 
    medianRaArray  = np.zeros(noOfGal)
    medianDecArray = np.zeros(noOfGal)

    ##################### updated by HJT#################
    medianResidualRaArray  = np.zeros(noOfGal)
    medianResidualDecArray = np.zeros(noOfGal) 
    #####################################################

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
    objIDs  = table.col('obj_id')[maskForGal]
    #print "len(objIDs)", len(objIDs)
    #print "objIDs.size", objIDs.size
    raObjs  = table.col('ra')[maskForGal]
    decObjs = table.col('dec')[maskForGal]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raObjs, decObjs, status = s2t.ds2tp(raObjs, decObjs, ra0, dec0) # transform the (RA, DEC)/degree into (xi, eta)/degree, but here (xi, eta) denoted by (raObjs, decObjs)
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    mjdObjs = table.col('mjd')[maskForGal]

    ##################### updated by HJT#################
    raErrObjs  = table.col('raErr')[maskForGal]  
    decErrObjs = table.col('decErr')[maskForGal] 
    #####################################################

    try:
        # store the first object values in here
        currObjID = objIDs[0]
        currRaObj  = [raObjs[0]]
        currDecObj = [decObjs[0]]
	##################### updated by HJT#################
    	currRaErrObj  = [raErrObjs[0]]
        currDecErrObj = [decErrObjs[0]] 
    	#####################################################
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
    	    ##################### updated by HJT#################
    	    raErr  = raErrObjs[i]  
    	    decErr = decErrObjs[i] 
    	    #####################################################

            if debug: print 'object ID',objID
            if debug: print 'ra of object',ra
            if debug: print 'dec of object',dec
            if debug: print 'mjd of object',mjd
            if (objID == currObjID):
                currRaObj.append(ra)
                currDecObj.append(dec)
                currMjdObj.append(mjd)
    	        ##################### updated by HJT#################
    	        currRaErrObj.append(raErr)
                currDecErrObj.append(decErr) 
    	        #####################################################

                if debug: print 'same object now'
                noOfObsArray[objIDcounter] += 1
            else:
                if debug: print 'different object now'
                if debug: print 'run functions for the', objIDcounter, 'object'
		currMjdSorted  = np.sort(currMjdObj)
    		deltaT = currMjdSorted[0:-1] - currMjdSorted[1:]
    		tBreakAt   = np.where(deltaT < (-100.0))[0]
    		currMjdBreakAt = currMjdSorted[tBreakAt+1]
		if debug: print "currMjdObj:", currMjdObj
                if(len(currMjdBreakAt) >= 3.0):#(noOfObsArray[objIDcounter] >= 2.0):
		    #print "tian test1,currMjdBreakAt", currMjdBreakAt
                    if debug: print 'no of obs ', noOfObsArray[objIDcounter]
                    # convert list into array
                    currRaObj = np.array(currRaObj)
                    currDecObj = np.array(currDecObj)
		    ##################### updated by HJT#################
		    currMjdObj = np.array(currMjdObj)
    	            currRaErrObj = np.array(currRaErrObj)
                    currDecErrObj = np.array(currDecErrObj)
		    '''currRaObj = currRaObj[currMjdObj<=currMjdBreakAt[2]]  ## Here notice we only choose the first 4 epochs 
		    currDecObj = currDecObj[currMjdObj<=currMjdBreakAt[2]]
		    currRaErrObj = currRaErrObj[currMjdObj<=currMjdBreakAt[2]]
		    currDecErrObj = currDecErrObj[currMjdObj<=currMjdBreakAt[2]]
		    currMjdObj = currMjdObj[currMjdObj<=currMjdBreakAt[2]]'''
    	            #####################################################

                    if debug: print 'current objID ra dec mjd :',currObjID, currRaObj, currDecObj, currMjdObj
                    # calculate and store medians
                    #medianRa = np.median(currRaObj)
                    #medianDec = np.median(currDecObj)
		    ##################### updated by HJT#################
    	            medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
                    medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	            #####################################################

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
		    indTmp = np.arange(pos,i)
		    #indTmp = indTmp[currMjdObj<=currMjdBreakAt[2]]
                    residualRaArray[indTmp] = residualRa
                    residualDecArray[indTmp] = residualDec
		    ##################### updated by HJT#################
		    # store median value for each object
                    medianResidualRaArray[objIDcounter] = np.median(residualRa)
                    medianResidualDecArray[objIDcounter] = np.median(residualDec)
		    #####################################################

                # store the objID no matter what
                objIDarray[objIDcounter] = currObjID
                objIDcounter += 1
                currObjID = objID
                pos = i
                # start new ra, dec, and mjd arrays and
                # don't forget to increment the noOfObsArray for this object
                currRaObj = [raObjs[i]]
                currDecObj = [decObjs[i]]
		##################### updated by HJT#################
    	        currRaErrObj = [raErrObjs[i]]
                currDecErrObj = [decErrObjs[i]]
    	        #####################################################
                currMjdObj = [mjdObjs[i]]
                noOfObsArray[objIDcounter] += 1

        # after going through all of the rows, make sure to process the last
        # object, if it has enough observations
        i +=1

        print 'processing the last object now'
        currMjdSorted  = np.sort(currMjdObj)
    	deltaT = currMjdSorted[0:-1] - currMjdSorted[1:]
    	tBreakAt   = np.where(deltaT < (-100.0))[0]
    	currMjdBreakAt = currMjdSorted[tBreakAt+1]
        if(len(currMjdBreakAt) >= 3.0):#(noOfObsArray[objIDcounter] >= 2.0):
            if debug: print 'obs greater than three, adding the last object too'
            currRaObj = np.array(currRaObj)
            currDecObj = np.array(currDecObj)
	    ##################### updated by HJT#################
    	    currRaErrObj = np.array(currRaErrObj)
            currDecErrObj = np.array(currDecErrObj)
	    currMjdObj = np.array(currMjdObj)
	    '''currRaObj = currRaObj[currMjdObj<=currMjdBreakAt[2]]  ## Here notice we only choose the first 4 epochs 
	    currDecObj = currDecObj[currMjdObj<=currMjdBreakAt[2]]
	    currRaErrObj = currRaErrObj[currMjdObj<=currMjdBreakAt[2]]
	    currDecErrObj = currDecErrObj[currMjdObj<=currMjdBreakAt[2]]
	    currMjdObj = currMjdObj[currMjdObj<=currMjdBreakAt[2]]'''
    	    #####################################################
            # calculate and store medians
            #medianRa = np.median(currRaObj)
            #medianDec = np.median(currDecObj)

	    ##################### updated by HJT#################
    	    medianRa = np.sum(currRaObj/currRaErrObj**2)/np.sum(1/currRaErrObj**2)
            medianDec = np.sum(currDecObj/currDecErrObj**2)/np.sum(1/currDecErrObj**2)
    	    #####################################################
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
            indTmp = np.arange(pos,i)
	    #indTmp = indTmp[currMjdObj<=currMjdBreakAt[2]]
            residualRaArray[indTmp] = residualRa
            residualDecArray[indTmp] = residualDec
	    ##################### updated by HJT#################
	    # store median value for each object
            medianResidualRaArray[objIDcounter] = np.median(residualRa)
            medianResidualDecArray[objIDcounter] = np.median(residualDec)
	    #####################################################
            # also need to store the objID, no matter what
            objIDarray[objIDcounter] = currObjID
       
        print 'total time taken is', time() - t1


    except ValueError:
        print "currObjID medianRa medianDec residualRa residualDec, pos, i, objIDcounter", currObjID, medianRa, medianDec, residualRa, residualDec, pos, i, objIDcounter
        pass
        
    return objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray, medianResidualRaArray, medianResidualDecArray




def calcMedianAndResiduals3(table, refDB, ra0, dec0, gnomonic, debug):

    #maskForGal  = table.get_where_list('sdss_type==3')
    maskForStar = table.get_where_list('gal==0')

    galIDs = table.col('obj_id')[maskForGal]
    raGal   = table.col('ra')[maskForGal]
    decGal  = table.col('dec')[maskForGal]
    galMjd = table.col('mjd')[maskForGal]

    galIDfinal, raFinalGal, decFinalGal, medianRaError, medianDecError = readDB(refDB)
    mjdSorted  = np.sort(galMjd)
    deltaT     = mjdSorted[0:-1] - mjdSorted[1:]
    tBreakAt   = np.where(deltaT < (-100.0))[0]
    mjdBreakAt = mjdSorted[tBreakAt + 1]
    #print mjdBreakAt

    distTmp = sphdist(ra0, dec0, raFinalGal, decFinalGal)
    angSepMask = (distTmp) <= (60.0/60.0)
    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
	#print "len(inde) in gal:", inde.size
	if(inde.size<10): print "Error is confusing"#return pixelNo, flag1, inde.size
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]

    except IndexError:
	print "Error is confusing"
        #return pixelNo, flag1, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    #galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	#print "ID", ID
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
        else:
                #make them numpy arrays
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]

    return offsetRaArray, offsetDecArray

#################### Parallelly for the all sky offset ##################
def calcMedianAndResiduals4(parameterListForPixel):

    pixelNo, pixRa, pixDec, galIDs, raGal, decGal, mmjd, \
	galIDfinal, raFinalGal, decFinalGal, medianRaError, medianDecError, gnomonic = parameterListForPixel

    distTmp = sphdist(pixRa, pixDec, raFinalGal, decFinalGal)
    angSepMask = (distTmp) <= (60.0/60.0)
    numInRadius = sum(angSepMask)
    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):# transform the (RA, DEC)/degree into (xi, eta)/degree
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, pixRa, pixDec)
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    #print "gal in the circle:", numInRadius

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
	#print "len(inde) in gal:", inde.size
	if(inde.size<10): return pixelNo, flag1, inde.size
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, pixRa, pixDec)
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]

    except IndexError:
        return pixelNo, flag1, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    #galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	#print "ID", ID
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
        else:
                #make them numpy arrays
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
		#print currentRaGal, raFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]

    offsetRaArray = offsetRaArray*1000*3600
    offsetDecArray = offsetDecArray*1000*3600

    moffsetRaInPixel = np.median(offsetRaArray)
    moffsetDecInPixel = np.median(offsetDecArray)

    eoffsetRaInPixel = (0.741*(np.percentile(offsetRaArray, 75) - np.percentile(offsetRaArray, 25)))
    eoffsetDecInPixel = (0.741*(np.percentile(offsetDecArray, 75) - np.percentile(offsetDecArray, 25)))
    #print pixelNo, pixRa, pixDec, moffsetRaInPixel, moffsetDecInPixel, eoffsetRaInPixel, eoffsetDecInPixel, numInRadius
    return pixelNo, pixRa, pixDec, mmjd, moffsetRaInPixel, moffsetDecInPixel, eoffsetRaInPixel, eoffsetDecInPixel, numInRadius

def offsetPixels(parameterListForPixel, datname):
    dat  = np.zeros(len(parameterListForPixel), dtype=[('pixelNo', 'i4'), ('ra', 'f8'), ('dec', 'f8'), ('mmjd', 'i4'), \
	('moffsetRaInPixel', 'f4'), ('moffsetDecInPixel', 'f4'), ('eoffsetRaInPixel', 'f4'), ('eoffsetDecInPixel', 'f4'), ('numInRadius', 'i4')])
    #start workers
    pool = Pool(processes=10)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(calcMedianAndResiduals4, parameterListForPixel[0:], chunksize=10)
    noIterated = 0
    for res in iterator:
        if (len(res)==3):
	    print "there was an error with pandas :( ", res
	else:
	    pixelNo_, ra_, dec_, mmjd_, moffsetRaInPixel_, moffsetDecInPixel_, eoffsetRaInPixel_, eoffsetDecInPixel_, numInPixel_ = res
	    dat[noIterated]['pixelNo'] = int(pixelNo_)
	    dat[noIterated]['ra'] = ra_
	    dat[noIterated]['dec'] = dec_
	    dat[noIterated]['mmjd'] = mmjd_
	    dat[noIterated]['moffsetRaInPixel'] = moffsetRaInPixel_
	    dat[noIterated]['moffsetDecInPixel'] = moffsetDecInPixel_
	    dat[noIterated]['eoffsetRaInPixel'] = eoffsetRaInPixel_
	    dat[noIterated]['eoffsetDecInPixel'] = eoffsetDecInPixel_
	    dat[noIterated]['numInRadius'] = int(numInPixel_)
	    noIterated += 1
    pool.terminate()
    dat = dat[0:noIterated]
    np.save(datname, dat)
    print "total pixel nums:", noIterated

##################################################################
def pixelTasks2(parameterListForPixel):
    """Workers will execute this function."""
    #unpack the input parameter list (i.e., separate them into different arrays)
    #objID -- contains the good galaxy IDs -- unique no of times
    #objIDs -- contains all the galaxy IDs corresponding to all the detections of galaxies
    #medRa and medDec -- contain good median values 
    #resRa and resDec -- contain the residual values for each detection of galaxy-- not good
    #ra/decAll -- ra/dec values for all the galaxies -- not good
    #pixAll -- pixel indices for all the galaxies -- not good
    pixelNo, pixRa, pixDec, objID, medRa, medDec, objIDs, MJDs, resRa, resDec, raAll, decAll, pixAll,mjdBreakAt, gnomonic, ra0, dec0 = parameterListForPixel


    #find objects within the searchRadius
    #angSepMask = sphdist(pixRa, pixDec, medRa, medDec) <= (searchRadius/60.0)
    distTmp = sphdist(pixRa, pixDec, medRa, medDec)
    angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>600):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:NGAL]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    objInRadius = objID[angSepMask]
    #print "galaxy number in the circle:",len(objInRadius) 
    #print "min(decAll), max(decAll):", min(decAll), max(decAll)

    '''we use PANDAS'''
    try:
        flag1 = 1 
	r1 = time()
	temp = pd.Index(objIDs)
	inde = temp.get_indexer_for(objInRadius)
	#print 'Test1:',max(inde), len(inde)
	inde = inde[inde > -1]
	#print 'Test2:',len(inde)
	inde = np.array(inde, dtype='i8')
	#print 'Test3:',resRa
	#print 'Test3:',len(resRa), max(inde)
	resRaValuesInRadius = resRa[inde]
	#print 'Test4:',resRaValuesInRadius
	resDecValuesInRadius = resDec[inde]
	mjdInRadius = MJDs[inde]
	#print 'galaxy mjdInRadius:',mjdInRadius
	    
    except IndexError:
	#print 'galaxy mjdInRadius:',mjdInRadius
        return pixelNo, flag1, inde.size
	    
    #for storing median values for each epoch
    medianResidualRaEpochWise  = np.zeros(len(mjdBreakAt)+1)
    medianResidualDecEpochWise = np.zeros(len(mjdBreakAt)+1)
    indexInPixel = pixAll == pixelNo
    objIDinPixel = objID[indexInPixel]
   
    try :
	    flag2 = 2 
	    temp2 = pd.Index(objIDs)
	    inde2 = temp2.get_indexer_for(objIDinPixel)
	    inde2 = inde2[inde2 > -1]
	    inde2 = np.array(inde2, dtype='i8')
	    raInPixel = raAll[inde2]
	    decInPixel = decAll[inde2]
	    #########################(RA, DEC) ==> (xi, eta)#################################
	    if(gnomonic):
	        raInPixel, decInPixel, status = s2t.ds2tp(raInPixel, decInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	        #print "Number points in bad transformation:", np.sum(status)
	    #################################################################################
	    objInPixel = objIDs[inde2]
	    mjdInPixel = MJDs[inde2]
	    #print "time taken for pandas indexer", time() - r1
    except IndexError:
	    return pixelNo, flag2, inde2.size
    
    #if there are atleast three epochs present
    if (len(mjdBreakAt)>=3.0):
	for mjdIdx in range(0, len(mjdBreakAt)+1):
	    if(mjdIdx==0):
	        var = mjdBreakAt[mjdIdx]
		mjdIndexInPixel = mjdInPixel < var
		mjdIndex = mjdInRadius < var
	    elif(mjdIdx==len(mjdBreakAt)):
	        var = mjdBreakAt[mjdIdx-1]
		mjdIndexInPixel = mjdInPixel >= var
		mjdIndex = mjdInRadius >= var
	    else:
	        var = mjdBreakAt[mjdIdx]
		mjdIndex = (mjdInRadius >= mjdBreakAt[mjdIdx-1]) & (mjdInRadius < var)
                mjdIndexInPixel = (mjdInPixel >= mjdBreakAt[mjdIdx-1]) & (mjdInPixel < var)

            if(any(mjdIndex)):
                resRaValues  = resRaValuesInRadius[mjdIndex]
                resDecValues = resDecValuesInRadius[mjdIndex]
		#print "number in resRaValues:", len(resRaValues)
                assert(resRaValues.size > 0.0), "resRaValues is empty"
                assert(resDecValues.size > 0.0), "resDecValues is empty"
                medianResidualRaEpochWise[mjdIdx]  = np.median(resRaValues)
                medianResidualDecEpochWise[mjdIdx] = np.median(resDecValues)
		#print "np.median(resRaValues):", np.median(resRaValues)
                raInPixel[mjdIndexInPixel] -=  medianResidualRaEpochWise[mjdIdx]
                decInPixel[mjdIndexInPixel] -=  medianResidualDecEpochWise[mjdIdx]

    #calculate final ra/dec as the median of newRa/newDec values
    finalRaArray  = np.zeros(objIDinPixel.size)
    finalDecArray = np.zeros(objIDinPixel.size)
    medianRaError  = np.zeros(objIDinPixel.size)
    medianDecError = np.zeros(objIDinPixel.size)
    
    objInPixel = np.array(objInPixel)
    raInPixel = np.array(raInPixel)
    decInPixel = np.array(decInPixel)
    mjdInPixel = np.array(mjdInPixel)
    #print "galaxy number in the pixel:", len(finalRaArray)
    try:
        #take the current values separately
        obj   = objInPixel[0]
        raObj = [raInPixel[0]]
        decObj = [decInPixel[0]]
        mjdObj = [mjdInPixel[0]]
        #set variables
        counter = 0
        pos2 = 0
        t2   = time()
        for i in np.arange(objInPixel.size-1)+1:#np.arange(10)+1: #
            #print 'on the row number',i
            objID = objInPixel[i]
            ra  = raInPixel[i]
            dec = decInPixel[i]
            if(objID == obj):
                raObj.append(ra)
                decObj.append(dec)
            else:
                raObj = np.array(raObj)
                decObj= np.array(decObj)

                finalRaArray[counter]  = np.median(raObj)
                finalDecArray[counter] = np.median(decObj)
                #calculate rms values for ra/dec
                rmsRa  = 0.741*(np.percentile(raObj, 75) - np.percentile(raObj, 25))
                rmsDec = 0.741*(np.percentile(decObj, 75) - np.percentile(decObj, 25))
                #calculate uncertainity in median coordinates
                medianRaError[counter]  = np.sqrt((np.pi/2)/(raObj.size-1))*rmsRa
                medianDecError[counter] = np.sqrt((np.pi/2)/(decObj.size-1))*rmsDec
                obj  = objID
                pos2 = i
                counter +=1
                #make them lists again
                raObj  =[raInPixel[i]]
                decObj =[decInPixel[i]]
        #processing the last object now
        raObj  = np.array(raObj)
        decObj = np.array(decObj)
        finalRaArray[counter]  = np.median(raObj)
        finalDecArray[counter] = np.median(decObj)
        #calculate rms values for ra/dec
        rmsRa  = 0.741*(np.percentile(raObj, 75) - np.percentile(raObj, 25))
        rmsDec = 0.741*(np.percentile(decObj, 75) - np.percentile(decObj, 25))
        medianRaError[counter]  = np.sqrt((np.pi/2)/(raObj.size-1))*rmsRa
        medianDecError[counter] = np.sqrt((np.pi/2)/(decObj.size-1))*rmsDec
	#########################(xi, eta) ==> (RA, DEC)#################################
	if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	    finalRaArray, finalDecArray = s2t.dtp2s(np.radians(finalRaArray), np.radians(finalDecArray), ra0, dec0)
	#################################################################################
        return objIDinPixel, finalRaArray, finalDecArray, medianRaError, medianDecError

    except IndexError:
        #pass
        print "index error"
        print pixelNo
        return np.array(objIDinPixel), np.array(finalRaArray), np.array(finalDecArray), np.array(medianRaError), np.array(medianDecError)


##################################################################
def referenceGals(parameterListForPixel, dbname):
    db = sqlite3.connect(dbname)
    cursor = db.cursor()
    cursor.execute('''DROP TABLE IF EXISTS finalData''')
    db.commit()  
    cursor = db.cursor()
    cursor.execute( '''CREATE TABLE finalData(objID INT NOT NULL, ra REAL NOT NULL, dec REAL NOT NULL, raError REAL NOT NULL, decError REAL NOT NULL)''')
    db.commit() 
    #start workers
    pool = Pool(processes=15)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(pixelTasks2, parameterListForPixel[0:], chunksize=300)
    counter = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        print(res)
        if (len(res)==3):
		print "there was an error with pandas :( ", res
	else:
		objID_, ra_, dec_, rErr_, dErr_ = res
		noOfPixelsIterated += 1
		print "objID in the pixel:",objID_
		#print "noOfPixelsIterated",noOfPixelsIterated
        #so that you donot encounter empty pixels! 
		if (objID_.size != 0):
			tmp = [(int(objID),float(ra), float(dec), float(rErr), float(dErr)) for objID, ra, dec, rErr, dErr in \
			    zip(objID_, ra_, dec_, rErr_, dErr_)]
			dat = dat + tmp
			counter = counter + len(objID_)
			'''if ((noOfPixelsIterated>50000)&(noOfPixelsIterated>int(noOfPixelsIterated/50000)*50000)): #
				cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', tmp)
				db.commit()
				#print "committing last" '''
			if (noOfPixelsIterated % 5000 == 0):
				cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
				db.commit()
				print "committing in 50000*N"
				dat = []
		elif (objID_.size ==0):
			print "zero array"   
		#print "counter", counter
    cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
    db.commit()
    print "counter=total no of rows added to database", counter
    #terminate the pool of multiprocessors
    pool.terminate()
    # close connection to db
    db.close()

def readDB(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size
	objIDarray = np.zeros(100000000, dtype = 'i8') # eight zeros
	finalRaArray  = np.zeros(100000000)
	finalDecArray = np.zeros(100000000)
	medianRaError = np.zeros(100000000)
	medianDecError = np.zeros(100000000)

	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT * FROM finalData'):
	    objIDarray[rowCounter]     = row[0]
	    finalRaArray[rowCounter]  = row[1]
	    finalDecArray[rowCounter]  = row[2]
	    medianRaError[rowCounter]  = row[3]
	    medianDecError[rowCounter]  = row[4]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	objIDarray = objIDarray[0:rowCounter]
	finalRaArray = finalRaArray[0:rowCounter]
	finalDecArray = finalDecArray[0:rowCounter]
	medianRaError = medianRaError[0:rowCounter]*1000*3600
	medianDecError = medianDecError[0:rowCounter]*1000*3600
        return objIDarray, finalRaArray, finalDecArray, medianRaError, medianDecError

def readPMDB(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size
	'''objID = np.zeros(100000000, dtype = 'i8')
	ra  = np.zeros(100000000)
	dec  = np.zeros(100000000)
	rmag  = np.zeros(100000000)
	muXI  = np.zeros(100000000)
	muErrXI  = np.zeros(100000000)
	muXxi  = np.zeros(100000000)
	muErrXxi  = np.zeros(100000000)
	muXIog  = np.zeros(100000000)
	muErrXIog  = np.zeros(100000000)
	chi2XI  = np.zeros(100000000)
	muETA  = np.zeros(100000000)
	muErrETA  = np.zeros(100000000)
	muXeta  = np.zeros(100000000)
	muErrXeta  = np.zeros(100000000)
	muETAog  = np.zeros(100000000)
	muErrETAog  = np.zeros(100000000)
	chi2ETA  = np.zeros(100000000)
	Nobs  = np.zeros(100000000)
	flag  = np.zeros(100000000)'''

	GPS1 = np.zeros(100000000, dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), 
            ('mr', 'f4'), ('muXI', 'f4'), ('muErrXI', 'f4'), ('muXxi', 'f4'), \
            ('muErrXxi', 'f4'), ('muXIog', 'f4'), ('muErrXIog', 'f4'), ('chi2XI', 'f4'), \
	    ('muETA', 'f4'), ('muErrETA', 'f4'), ('muXeta', 'f4'), ('muErrXeta', 'f4'), \
	    ('muETAog', 'f4'), ('muErrETAog', 'f4'), ('chi2ETA', 'f4'),('Nobs', 'u2'), ('flag', 'u1')]) 

	names = GPS1.dtype.names
	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT * FROM finalData'):
	    '''objID[rowCounter]     = row[0]
	    ra[rowCounter]  = row[1]
	    dec[rowCounter]  = row[2]
	    rmag[rowCounter]  = row[3]
	    muXI[rowCounter]  = row[4]
	    muErrXI[rowCounter]  = row[5]
	    muXxi[rowCounter]  = row[6]
	    muErrXxi[rowCounter]  = row[7]
	    muXIog[rowCounter]  = row[8]
	    muErrXIog[rowCounter]  = row[9]
	    chi2XI[rowCounter]  = row[10]
	    muETA[rowCounter]  = row[11]
	    muErrETA[rowCounter]  = row[12]
	    muXeta[rowCounter]  = row[13]
	    muErrXeta[rowCounter]  = row[14]
	    muETAog[rowCounter]  = row[15]
	    muErrETAog[rowCounter]  = row[16]
	    chi2ETA[rowCounter]  = row[17]
	    Nobs[rowCounter]  = row[18]
	    flag[rowCounter]  = row[19]'''
	    for idx in range(0, len(names)):
		GPS1[names[idx]][rowCounter] = row[idx]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	GPS1 = GPS1[0:rowCounter]
        return GPS1

def readPMDB2(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size


	GPS1 = np.zeros(100000000, dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), 
            ('mr', 'f4'), ('muXI', 'f4'), ('muErrXI', 'f4'), ('muXxi', 'f4'), ('muErrXxi', 'f4'), ('muXIog', 'f4'), \
	    ('muErrXIog', 'f4'), ('chi2XI', 'f4'), ('muETA', 'f4'), ('muErrETA', 'f4'), ('muXeta', 'f4'), \
	    ('muErrXeta', 'f4'), ('muETAog', 'f4'), ('muErrETAog', 'f4'), ('chi2ETA', 'f4'), ('Nobs', 'u2'), ('flag', 'u2')]) 

	names = GPS1.dtype.names
	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT objID, ra, dec, mr, muXI, muErrXI, muXxi, muErrXxi, muXIog, muErrXIog, \
	    chi2XI, muETA, muErrETA, muXeta, muErrXeta, muETAog, muErrETAog, chi2ETA, Nobs, flag FROM finalData'):
	    for idx in range(0, len(names)):
		GPS1[names[idx]][rowCounter] = row[idx]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	GPS1 = GPS1[0:rowCounter]
        return GPS1

def readPMDB3(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size


	GPS1 = np.zeros(100000000, dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), 
            ('mr', 'f4'), ('muXI', 'f4'), ('muErrXI', 'f4'), ('muXxi', 'f4'), ('muErrXxi', 'f4'), ('muXIog', 'f4'), \
	    ('muErrXIog', 'f4'), ('chi2XI', 'f4'), ('muETA', 'f4'), ('muErrETA', 'f4'), ('muXeta', 'f4'), \
	    ('muErrXeta', 'f4'), ('muETAog', 'f4'), ('muErrETAog', 'f4'), ('chi2ETA', 'f4'), \
	    ('muXIps1','f4'), ('muErrXIps1','f4'), ('chi2XIps1','f4'), ('muETAps1','f4'), ('muErrETAps1','f4'), \
	    ('chi2ETAps1','f4'), ('Nobs0','u2'), ('Nobs', 'u2'), ('flag', 'u2')]) 

	names = GPS1.dtype.names
	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT objID, ra, dec, mr, muXI, muErrXI, muXxi, muErrXxi, muXIog, muErrXIog, \
	    chi2XI, muETA, muErrETA, muXeta, muErrXeta, muETAog, muErrETAog, chi2ETA, muXIps1, muErrXIps1, chi2XIps1, \
	    muETAps1, muErrETAps1, chi2ETAps1, Nobs0, Nobs, flag FROM finalData'):
	    for idx in range(0, len(names)):
		GPS1[names[idx]][rowCounter] = row[idx]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	GPS1 = GPS1[0:rowCounter]
        return GPS1

def readPMDB4(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size


	GPS1 = np.zeros(100000000, dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), 
            ('mr', 'f4'), ('muXI', 'f4'), ('muErrXI', 'f4'), ('chi2XI', 'f4'), ('gpmra', 'f4'), ('gpmraErr', 'f4'), \
            ('muETA', 'f4'), ('muErrETA', 'f4'), ('chi2ETA', 'f4'), ('gpmdec', 'f4'), ('gpmdecErr', 'f4'), \
            ('Nobs0','u2'), ('Nobs', 'u2'), ('flag', 'u2')]) 

	names = GPS1.dtype.names
	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT objID, ra, dec, mr, muXI, muErrXI, \
	    chi2XI, gpmra, gpmraErr, muETA, muErrETA, chi2ETA, \
	    gpmdec, gpmdecErr, Nobs0, Nobs, flag FROM finalData'):
	    for idx in range(0, len(names)):
		GPS1[names[idx]][rowCounter] = row[idx]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	GPS1 = GPS1[0:rowCounter]
        return GPS1

def readPMDB5(dbname):
	#establish connection
	connection = sqlite3.connect(dbname)
	#obtain cursor object
	cursor = connection.cursor()
	#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
	rowCounter = 0 
	#define arrays with a reasonable size


	GPS1 = np.zeros(100000000, dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), 
            ('mr', 'f4'), ('muXI', 'f4'), ('muErrXI', 'f4'), ('muXxi', 'f4'), ('muErrXxi', 'f4'), ('chi2XI', 'f4'), \
            ('muXI_mcmc', 'f4'), ('muErrXI_mcmc', 'f4'), ('chi2XI_mcmc', 'f4'), ('gpmra', 'f4'), ('gpmraErr', 'f4'), \
            ('muETA', 'f4'), ('muErrETA', 'f4'), ('muXeta', 'f4'), ('muErrXeta', 'f4'), ('chi2ETA', 'f4'), \
            ('muETA_mcmc', 'f4'), ('muErrETA_mcmc', 'f4'), ('chi2ETA_mcmc', 'f4'), ('gpmdec', 'f4'), ('gpmdecErr', 'f4'), \
            ('Nobs0','u2'), ('Nobs', 'u2'), ('flag', 'u2')]) 

	names = GPS1.dtype.names
	#loop over all rows in the table and obtain values
	for row in cursor.execute('SELECT objID, ra, dec, mr, muXI, muErrXI, muXxi, muErrXxi, chi2XI, muXI_mcmc, muErrXI_mcmc, \
	    chi2XI_mcmc, gpmra, gpmraErr, muETA, muErrETA, muXeta, muErrXeta, chi2ETA, muETA_mcmc, muErrETA_mcmc, chi2ETA_mcmc, \
	    gpmdec, gpmdecErr, Nobs0, Nobs, flag FROM finalData'):
	    for idx in range(0, len(names)):
		GPS1[names[idx]][rowCounter] = row[idx]
	    rowCounter +=1

	#print "rowCounter", rowCounter
	#close the connection to database
	connection.close()

	#trim the array to the required size 
	GPS1 = GPS1[0:rowCounter]
        return GPS1

def npy2Table(npy, Star, h5file):
    h5file = tables.open_file(h5file, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    counterForRows = 0 
    dtype = table.colnames
    print dtype
    for idx in xrange(0, len(npy['obj_id'])):
	for j in range(len(dtype)): 
            star[dtype[j]] = npy[dtype[j]][idx]
        star.append()
	counterForRows = counterForRows + 1
        if (counterForRows % 50000 == 0): 
            table.flush()
	    #star['obj_id'][0:10]
    table.flush()
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()


def pixelTasksCombinedData(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd , starIDs, starMjds, raStar, decStar, rMagStar, nObsStar, pixAllStar, pixelIndexForStar_l,mjdBreakAt, gnomonic, ra0, dec0 = parameterList

    distTmp = sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>600):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:NGAL]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]

    except IndexError:
        return pixelNo, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    print "galIDinRadius:", galIDinRadius
    
    galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
                #print "same object now"
        else:
                #make them numpy arrays
		#print "currentGalID, uniqueGalIDinRadius:", currentGalID, uniqueGalIDinRadius[galCounter]
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
		#print currentRaGal, raFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
		#print "offsetRa, offsetDec", offsetRa, offsetDec
                #update counters and variables
	        galCounter+= 1
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]
    #print "time taken to calculate offsets of galaxies from their final ra/dec", time()-t
    print "galCounter after loop", galCounter	

	
    # select stars inside the pixel
   
    #resIndexInPixel =  pixAllStar == pixelNo
    #print "pixelNo:", pixelNo, pixAllStar
    #print "sum of star number in Pixel:", np.sum(resIndexInPixel)
    
    #########################################updata by Nie#########################
    resIndexInPixel= []
    for j in pixelIndexForStar_l:
        resIndexInPixel_l =  int(pixelNo*10**(-9)) == j
        resIndexInPixel.append(resIndexInPixel_l)
    #print pickPixelNo, pixelIndexForStar
        #if resIndexInPixel != False:
    print "star in Pixel:", np.sum(resIndexInPixel)
        #resIndexInPixel_all.append(resIndexInPixel)
    ########################################undata by Nie#########################
    
    #uniqueStarIDinPixel = starID[resIndexInPixel]
    #match all detections of stars with those in the pixel
    #indexInPixel = np.in1d(starIDs, uniqueStarIDinPixel)
    #obtain original ra and dec values for stars within pixel
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	#print "raStarInPixel, decStarInPixel:", raStarInPixel, decStarInPixel
	raStarInPixel, decStarInPixel, status = s2t.ds2tp(raStarInPixel, decStarInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    #raErrStarInPixel  = raErrStar[resIndexInPixel]
    #decErrStarInPixel = decErrStar[resIndexInPixel]
    rMagStarInPixel = rMagStar[resIndexInPixel]
    nObsStarInPixel = nObsStar[resIndexInPixel]
    #raStarErrInPixel  = np.zeros(np.sum(resIndexInPixel))
    #decStarErrInPixel = np.zeros(np.sum(resIndexInPixel))
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]
    print 'sum number of resIndexInPixel', np.sum(resIndexInPixel)

    medianRaOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    medianDecOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    #avgMjd = np.zeros(len(mjdBreakAt)+1)

    if (len(mjdBreakAt)>= 1.0):
        for mjdIdx in range(0, len(mjdBreakAt)+1):
	    if(mjdIdx==0):
	        var = mjdBreakAt[mjdIdx]
		mjdIndexInPixel = mjdStarInPixel < var
		mjdIndex = galMjdInRadius < var
	    elif(mjdIdx==len(mjdBreakAt)):
	        var = mjdBreakAt[mjdIdx-1]
		mjdIndexInPixel = mjdStarInPixel >= var
		mjdIndex = galMjdInRadius >= var
	    else:
	        var = mjdBreakAt[mjdIdx]
		mjdIndex = (galMjdInRadius >= mjdBreakAt[mjdIdx-1]) & (galMjdInRadius < var)
                mjdIndexInPixel = (mjdStarInPixel >= mjdBreakAt[mjdIdx-1]) & (mjdStarInPixel < var)
            if (any(mjdIndex)):
                offsetRaValues  = offsetRaArray[mjdIndex]
                offsetDecValues = offsetDecArray[mjdIndex]
                medianRaOffsetEpochWise[mjdIdx]  = np.median(offsetRaValues)
                medianDecOffsetEpochWise[mjdIdx] = np.median(offsetDecValues)
                #update the ra/dec - replace old values
                raStarInPixel[mjdIndexInPixel] -= medianRaOffsetEpochWise[mjdIdx]
                decStarInPixel[mjdIndexInPixel]-= medianDecOffsetEpochWise[mjdIdx]
		#print "np.median(offsetRaValues):", np.median(offsetRaValues)
		if((offsetRaValues.size>1)):
		    raErrTmp = np.sqrt((np.pi/2)/(offsetRaValues.size-1)) * \
		        (0.741*(np.percentile(offsetRaValues, 75) - np.percentile(offsetRaValues, 25)))
		    decErrTmp = np.sqrt((np.pi/2)/(offsetDecValues.size-1)) * \
		        (0.741*(np.percentile(offsetDecValues, 75) - np.percentile(offsetDecValues, 25)))
		    #raErrStarInPixel[mjdIndexInPixel] = np.sqrt(raErrStarInPixel[mjdIndexInPixel]**2 + raErrTmp**2) # 
		    #decErrStarInPixel[mjdIndexInPixel] = np.sqrt(decErrStarInPixel[mjdIndexInPixel]**2 + decErrTmp**2) # 
                #avgMjd[var1] = (mjdStarInPixel[mjdIndexInPixel].max() - mjdStarInPixel[mjdIndexInPixel].min() )/2

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel

    #########################(xi, eta) ==> (RA, DEC)#################################
    if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(finalRaArrayStar), np.radians(finalDecArrayStar), ra0, dec0)
    #################################################################################
    if(len(finalObjIDstar)<1):
	return pixelNo, inde.size
    else:
        return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, mjdStarInPixel,rMagStarInPixel, nObsStarInPixel


def pixelTasksCombinedDataSDSS(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd , starIDs, starMjds, raStar, decStar, raErrStar, decErrStar, pixAllStar, mjdBreakAt, gnomonic, ra0, dec0 = parameterList
    distTmp = sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>200):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:200]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]
	
    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
	#print "len(inde) in gal:", inde.size
	if(inde.size<10): return pixelNo, flag1, inde.size
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]
    except IndexError:
        return pixelNo, flag1, inde.size
      
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",len(uniqueGalIDinRadius), noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
                #print "same object now"
        else:
                #make them numpy arrays
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
		#print currentRaGal, raFinalGalInRadius[uniqueGalIDinRadius==currentGalID], decFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
		#print "offsetRa, offsetDec", offsetRa, offsetDec
                #update counters and variables
	        galCounter+= 1
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]
    #print "time taken to calculate offsets of galaxies from their final ra/dec", time()-t
    #print "galCounter after loop", galCounter	

	
    # select stars inside the pixel
    resIndexInPixel =  pixAllStar == pixelNo
    #print "pixelNo:", pixelNo
    #print "sum of star number in Pixel:", np.sum(resIndexInPixel)
    #uniqueStarIDinPixel = starID[resIndexInPixel]
    #match all detections of stars with those in the pixel
    #indexInPixel = np.in1d(starIDs, uniqueStarIDinPixel)
    #obtain original ra and dec values for stars within pixel
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	#print "raStarInPixel, decStarInPixel:", raStarInPixel, decStarInPixel
	raStarInPixel, decStarInPixel, status = s2t.ds2tp(raStarInPixel, decStarInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    raErrStarInPixel  = raErrStar[resIndexInPixel]
    decErrStarInPixel = decErrStar[resIndexInPixel]
    #rMagStarInPixel = rMagStar[resIndexInPixel]
    #raStarErrInPixel  = np.zeros(np.sum(resIndexInPixel))
    #decStarErrInPixel = np.zeros(np.sum(resIndexInPixel))
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]
    #print 'sum number of resIndexInPixel', np.sum(resIndexInPixel)

    medianRaOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    medianDecOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    #avgMjd = np.zeros(len(mjdBreakAt)+1)

    if (len(mjdBreakAt)>= 1.0):
        for mjdIdx in range(0, len(mjdBreakAt)+1):
	    if(mjdIdx==0):
	        var = mjdBreakAt[mjdIdx]
		mjdIndexInPixel = mjdStarInPixel < var
		mjdIndex = galMjdInRadius < var
	    elif(mjdIdx==len(mjdBreakAt)):
	        var = mjdBreakAt[mjdIdx-1]
		mjdIndexInPixel = mjdStarInPixel >= var
		mjdIndex = galMjdInRadius >= var
	    else:
	        var = mjdBreakAt[mjdIdx]
		mjdIndex = (galMjdInRadius >= mjdBreakAt[mjdIdx-1]) & (galMjdInRadius < var)
                mjdIndexInPixel = (mjdStarInPixel >= mjdBreakAt[mjdIdx-1]) & (mjdStarInPixel < var)
            if (any(mjdIndex)):
                offsetRaValues  = offsetRaArray[mjdIndex]
                offsetDecValues = offsetDecArray[mjdIndex]
                medianRaOffsetEpochWise[mjdIdx]  = np.median(offsetRaValues)
                medianDecOffsetEpochWise[mjdIdx] = np.median(offsetDecValues)
                #update the ra/dec - replace old values
                raStarInPixel[mjdIndexInPixel] -= medianRaOffsetEpochWise[mjdIdx]
                decStarInPixel[mjdIndexInPixel]-= medianDecOffsetEpochWise[mjdIdx]
		#print "np.median(offsetRaValues):", np.median(offsetRaValues)
		if((offsetRaValues.size>1)):
		    raErrTmp = np.sqrt((np.pi/2)/(offsetRaValues.size-1)) * \
		        (0.741*(np.percentile(offsetRaValues, 75) - np.percentile(offsetRaValues, 25)))
		    decErrTmp = np.sqrt((np.pi/2)/(offsetDecValues.size-1)) * \
		        (0.741*(np.percentile(offsetDecValues, 75) - np.percentile(offsetDecValues, 25)))
		    #print "raErrStarInPixel[mjdIndexInPixel], offsetRaValues:", raErrStarInPixel[mjdIndexInPixel], offsetRaValues
		    raErrStarInPixel[mjdIndexInPixel] = np.sqrt(raErrStarInPixel[mjdIndexInPixel]**2 + raErrTmp**2)
		    decErrStarInPixel[mjdIndexInPixel] = np.sqrt(decErrStarInPixel[mjdIndexInPixel]**2 + decErrTmp**2)
                #avgMjd[var1] = (mjdStarInPixel[mjdIndexInPixel].max() - mjdStarInPixel[mjdIndexInPixel].min() )/2

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel

    #########################(xi, eta) ==> (RA, DEC)#################################
    if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(finalRaArrayStar), np.radians(finalDecArrayStar), ra0, dec0)
    #################################################################################
    #print "len(finalObjIDstar):", len(finalObjIDstar)
    if(len(finalObjIDstar)<1):
	return pixelNo, flag1, inde.size
    else:
        return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel

def pixelTasksCombinedDataSDSS2(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd, starIDs, starMjds, raStar, decStar, raErrStar, decErrStar, pixAllStar, mjdBreakAt, gnomonic, ra0, dec0 = parameterList

    distTmp = sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>200):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:200]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
	#print "len(inde) in gal:", inde.size
	if(inde.size<10): return pixelNo, flag1, inde.size
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]

    except IndexError:
        return pixelNo, flag1, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    #galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",len(uniqueGalIDinRadius), noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	#print "ID", ID
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
                #print "same object now"
        else:
                #make them numpy arrays
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
		#print currentRaGal, raFinalGalInRadius[uniqueGalIDinRadius==currentGalID], decFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
		#print "offsetRa:", offsetRa
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
                #update counters and variables
	        #galCounter+= 1
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]
    #print "time taken to calculate offsets of galaxies from their final ra/dec", time()-t
    #print "galCounter after loop", galCounter	

	
    # select stars inside the pixel
    resIndexInPixel =  pixAllStar == pixelNo
    #print "pixelNo:", pixelNo
    #print "sum of star number in Pixel:", np.sum(resIndexInPixel)
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	#print "raStarInPixel, decStarInPixel:", raStarInPixel, decStarInPixel
	raStarInPixel, decStarInPixel, status = s2t.ds2tp(raStarInPixel, decStarInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    raErrStarInPixel  = raErrStar[resIndexInPixel]
    decErrStarInPixel = decErrStar[resIndexInPixel]
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]
    #print 'sum number of resIndexInPixel', np.sum(resIndexInPixel)


    medianRaOffsetEpochWise  = np.median(offsetRaArray)
    medianDecOffsetEpochWise = np.median(offsetDecArray)
    #update the ra/dec - replace old values
    raStarInPixel -= medianRaOffsetEpochWise
    decStarInPixel-= medianDecOffsetEpochWise
    if((offsetRaArray.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetRaArray.size-1)) * \
	    (0.741*(np.percentile(offsetRaArray, 75) - np.percentile(offsetRaArray, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetDecArray.size-1)) * \
	    (0.741*(np.percentile(offsetDecArray, 75) - np.percentile(offsetDecArray, 25)))
	#print "raErrStarInPixel, offsetRaArray:", raErrStarInPixel, offsetRaArray
	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel

    #########################(xi, eta) ==> (RA, DEC)#################################
    if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(finalRaArrayStar), np.radians(finalDecArrayStar), ra0, dec0)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo, flag1, inde.size
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel

def pixelTasksCombinedDataGAIADRT(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd, starIDs, starMjds, raStar, \
        decStar, raErrStar, decErrStar,pmraStar,pmdecStar,epmraStar,epmdecStar,pixAllStar, mjdBreakAt, gnomonic, ra0, dec0 = parameterList

    distTmp = sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>200):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:200]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
	if(inde.size<10): return pixelNo, flag1, inde.size
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]

    except IndexError:
        return pixelNo, flag1, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
        else:
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID][0]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]

    # select stars inside the pixel
    resIndexInPixel =  pixAllStar == pixelNo
    #print "pixelNo:", pixelNo
    #print "sum of star number in Pixel:", np.sum(resIndexInPixel)
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	#print "raStarInPixel, decStarInPixel:", raStarInPixel, decStarInPixel
	raStarInPixel, decStarInPixel, status = s2t.ds2tp(raStarInPixel, decStarInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    raErrStarInPixel  = raErrStar[resIndexInPixel]
    decErrStarInPixel = decErrStar[resIndexInPixel]
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]
    pmraStarInPixel = pmraStar[resIndexInPixel]
    pmdecStarInPixel = pmdecStar[resIndexInPixel]
    epmraStarInPixel = epmraStar[resIndexInPixel]
    epmdecStarInPixel = epmdecStar[resIndexInPixel]

    #print 'sum number of resIndexInPixel', np.sum(resIndexInPixel)


    medianRaOffsetEpochWise  = np.median(offsetRaArray)
    medianDecOffsetEpochWise = np.median(offsetDecArray)
    #update the ra/dec - replace old values
    raStarInPixel -= medianRaOffsetEpochWise
    decStarInPixel-= medianDecOffsetEpochWise
    if((offsetRaArray.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetRaArray.size-1)) * \
	    (0.741*(np.percentile(offsetRaArray, 75) - np.percentile(offsetRaArray, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetDecArray.size-1)) * \
	    (0.741*(np.percentile(offsetDecArray, 75) - np.percentile(offsetDecArray, 25)))
	#print "raErrStarInPixel, offsetRaArray:", raErrStarInPixel, offsetRaArray
	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel

    #########################(xi, eta) ==> (RA, DEC)#################################
    if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(finalRaArrayStar), np.radians(finalDecArrayStar), ra0, dec0)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo, flag1, inde.size
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel,pmraStarInPixel,pmdecStarInPixel,epmraStarInPixel,epmdecStarInPixel



def pixelTasksCombinedDataGAIA(parameterList):
    # using X-validation, so the speed is too slow
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, pixelIndexForObj, \
	CRA, CDEC = parameterList

    distTmp = sphdist(pixelRa, pixelDec, graRef, gdecRef)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>600):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:NGAL]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    cN = len(angSepMask)
    objInRadius = gidRef[angSepMask]


    xi0, eta0, status = s2t.ds2tp(ra0, dec0, CRA, CDEC)
    xi1, eta1, status = s2t.ds2tp(ra1, dec1, CRA, CDEC)
    xi2, eta2, status = s2t.ds2tp(ra2, dec2, CRA, CDEC)
    gxi, geta, status = s2t.ds2tp(gra, gdec, CRA, CDEC)

    try:
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	for idx1 in range(0, cN):
	    ind0 = obj_id0 == objInRadius[idx1]
	    xtmp0 = mjd0[ind0]
	    xitmp0 = xi0[ind0]
	    etatmp0 = eta0[ind0]

	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1 = mjd1[ind1]
	    xitmp1 = xi1[ind1]
	    etatmp1 = eta1[ind1]	

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2 = mjd2[ind2]
	    xitmp2 = xi2[ind2]
	    etatmp2 = eta2[ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)

	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	        raErr2[ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	        decErr2[ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2)):
	        num0 = len(raErr0[ind0]) # Here num0 means the num of PS1		
	        xpm = np.zeros([8, num0])
	        for idxX in range(0, num0): # X means cross-validation
		    yErrtmpXI[idxX] = yErrtmpXI[idxX] + 9999.0
		    yErrtmpETA[idxX] = yErrtmpETA[idxX] + 9999.0
		    pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		    pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		    xpm[0, idxX] = pX[0]*365.25
		    xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
		    xpm[2, idxX] = abs(ytmpXI[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		    xpm[3, idxX] = pX[1]
		    xpm[4, idxX] = pXeta[0]*365.25
		    xpm[5, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
		    xpm[6, idxX] = abs(ytmpETA[idxX] - (pXeta[0]*xtmp[idxX] + pXeta[1]))
		    xpm[7, idxX] = pXeta[1]
		    yErrtmpXI[idxX] = yErrtmpXI[idxX] - 9999.0
		    yErrtmpETA[idxX] = yErrtmpETA[idxX] - 9999.0

		sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
	        idxO = abs(xpm[2] - np.median(xpm[2]))>3.0*sigTmp
		sigTmp_ = (0.741*(np.percentile(xpm[5], 75) - np.percentile(xpm[5], 25)))
	        idxO_ = abs(xpm[5] - np.median(xpm[5]))>3.0*sigTmp_
	        #print "sum(idxO), sum(idxO_):", sum(idxO), sum(idxO_)

    
	        if((sum(idxO)) & (len(xtmp)>5)):
		    idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
		    muXxi = xpm[0][idxO]
	            slpxi = xpm[3][idxO]
	        else:
		    p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		    muXxi = p[0]*365.0
		    slpxi = p[1]


	        if((sum(idxO_)) & (len(xtmp)>5)):
		    idxO_ = np.argmax(abs(xpm[6]))  # the max delta are probably a outlier
		    muXeta = xpm[4][idxO_]
	            slpeta = xpm[7][idxO_]
	        else:
		    peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		    muXeta = peta[0]*365.0
		    slpeta = peta[1]

	        gind = gobj_id == objInRadius[idx1]
	        gxtmp = gmjd[gind]
	        gxitmp = gxi[gind]
	        getatmp = geta[gind]

	        offsetXi[idx1] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	        offsetEta[idx1] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
		#print idx1, offsetXi[idx1], offsetEta[idx1]

    except IndexError:
        return pixelNo


    indInPixel = pixelIndexForObj == pixelNo
    print "sum(indInPixel):", pixelNo, sum(indInPixel)

    xiStarInPixel  = gxi[indInPixel]
    etaStarInPixel = geta[indInPixel]

    raErrStarInPixel  = graErr[indInPixel]
    decErrStarInPixel = gdecErr[indInPixel]
    starIDInPixel  = gobj_id[indInPixel]
    mjdStarInPixel = gmjd[indInPixel]

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel


def pixelTasksCombinedDataGAIA2(parameterList):
    # Do not use X-validation, to improve the speed
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, pixelIndexForObj, \
	CRA, CDEC = parameterList

    distTmp = sphdist(pixelRa, pixelDec, graRef, gdecRef)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>600):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:600]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]

    cN = 100#len(angSepMask)
    objInRadius = gidRef[angSepMask]


    xi0, eta0, status = s2t.ds2tp(ra0, dec0, CRA, CDEC)
    xi1, eta1, status = s2t.ds2tp(ra1, dec1, CRA, CDEC)
    xi2, eta2, status = s2t.ds2tp(ra2, dec2, CRA, CDEC)
    gxi, geta, status = s2t.ds2tp(gra, gdec, CRA, CDEC)

    try:
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	for idx1 in range(0, cN):
	    ind0 = obj_id0 == objInRadius[idx1]
	    xtmp0 = mjd0[ind0]
	    xitmp0 = xi0[ind0]
	    etatmp0 = eta0[ind0]

	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1 = mjd1[ind1]
	    xitmp1 = xi1[ind1]
	    etatmp1 = eta1[ind1]	

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2 = mjd2[ind2]
	    xitmp2 = xi2[ind2]
	    etatmp2 = eta2[ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)

	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	        raErr2[ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	        decErr2[ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2)):	        
		p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXxi = p[0]*365.0
		slpxi = p[1]

		peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muXeta = peta[0]*365.0
		slpeta = peta[1]

	        gind = gobj_id == objInRadius[idx1]
	        gxtmp = gmjd[gind]
	        gxitmp = gxi[gind]
	        getatmp = geta[gind]

	        offsetXi[idx1] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	        offsetEta[idx1] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
		#print idx1, offsetXi[idx1], offsetEta[idx1]

    except IndexError:
        return pixelNo


    indInPixel = pixelIndexForObj == pixelNo
    print "sum(indInPixel):", pixelNo, sum(indInPixel)

    xiStarInPixel  = gxi[indInPixel]
    etaStarInPixel = geta[indInPixel]

    raErrStarInPixel  = graErr[indInPixel]
    decErrStarInPixel = gdecErr[indInPixel]
    starIDInPixel  = gobj_id[indInPixel]
    mjdStarInPixel = gmjd[indInPixel]

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel

def pixelTasksCombinedDataGAIA3(parameterList):
    # Do not use X-validation, to improve the speed
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, pixelIndexForObj, \
	CRA, CDEC = parameterList


    cellsize = 0.2 #deg
    tmpMask0 = (abs(ra0-pixelRa)<cellsize)&(abs(dec0-pixelDec)<cellsize)
    tmpMask1 = (abs(ra1-pixelRa)<cellsize)&(abs(dec1-pixelDec)<cellsize)
    tmpMask2 = (abs(ra2-pixelRa)<cellsize)&(abs(dec2-pixelDec)<cellsize)
    tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)
    tmpMask4 = (abs(gra-pixelRa)<cellsize)&(abs(gdec-pixelDec)<cellsize)
    #print sum(tmpMask0), sum(tmpMask1), sum(tmpMask2), sum(tmpMask3), sum(tmpMask4)

    Nref = sum(tmpMask3)
    print "Nref:", Nref

    distTmp = sphdist(pixelRa, pixelDec, graRef[tmpMask3], gdecRef[tmpMask3])
    angSepMask = np.argsort(distTmp)
    if(Nref>100): 
	Nref = 100
    angSepMask = angSepMask[0:Nref]
    cN = Nref#len(angSepMask)
    objInRadius = gidRef[tmpMask3][angSepMask]


    try:
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	for idx1 in range(0, cN):
	    ind0 = obj_id0[tmpMask0] == objInRadius[idx1]
	    xtmp0 = mjd0[tmpMask0][ind0]
	    xitmp0 = ra0[tmpMask0][ind0]
	    etatmp0 = dec0[tmpMask0][ind0]

	    ind1 = obj_id1[tmpMask1] == objInRadius[idx1]
	    xtmp1 = mjd1[tmpMask1][ind1]
	    xitmp1 = ra1[tmpMask1][ind1]
	    etatmp1 = dec1[tmpMask1][ind1]	

	    ind2 = obj_id2[tmpMask2] == objInRadius[idx1]
	    xtmp2 = mjd2[tmpMask2][ind2]
	    xitmp2 = ra2[tmpMask2][ind2]
	    etatmp2 = dec2[tmpMask2][ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)

	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[tmpMask0][ind0], raErr1[tmpMask1][ind1], \
	        raErr2[tmpMask2][ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[tmpMask0][ind0], decErr1[tmpMask1][ind1], \
	        decErr2[tmpMask2][ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2)):	        
		p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXxi = p[0]*365.0
		slpxi = p[1]

		peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muXeta = peta[0]*365.0
		slpeta = peta[1]

	        gind = gobj_id[tmpMask4] == objInRadius[idx1]
	        gxtmp = gmjd[tmpMask4][gind]
	        gxitmp = gra[tmpMask4][gind]
	        getatmp = gdec[tmpMask4][gind]

	        offsetXi[idx1] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	        offsetEta[idx1] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
		#print idx1, offsetXi[idx1], offsetEta[idx1]

    except IndexError:
        return pixelNo


    indInPixel = pixelIndexForObj[tmpMask4] == pixelNo
    print "sum(indInPixel):", pixelNo, sum(indInPixel)

    xiStarInPixel  = gra[tmpMask4][indInPixel]
    etaStarInPixel = gdec[tmpMask4][indInPixel]

    raErrStarInPixel  = graErr[tmpMask4][indInPixel]
    decErrStarInPixel = gdecErr[tmpMask4][indInPixel]
    starIDInPixel  = gobj_id[tmpMask4][indInPixel]
    mjdStarInPixel = gmjd[tmpMask4][indInPixel]

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel

def pixelTasksCombinedDataGAIA4(parameterList):
    # Do not use X-validation, to improve the speed
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, pixelIndexForObj, \
	CRA, CDEC = parameterList


    t0= time()

    cellsize = 0.5 #deg

    #tmpMask3 = sphdist(pixelRa, pixelDec, graRef, gdecRef)<cellsize
    tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize) 
    Ntmp = sum(tmpMask3)
    if((Ntmp<300)&(Ntmp>50)): 
	cellsize = 1.0
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)
	Ntmp = sum(tmpMask3)
    if((Ntmp<50)): 
	cellsize = 1.5
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)
    '''elif(Ntmp>5000): 
	cellsize = 0.1
        #tmpMask3 = sphdist(pixelRa, pixelDec, graRef, gdecRef)<cellsize
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)'''
    '''tmpMask0 = sphdist(pixelRa, pixelDec, ra0, dec0)<cellsize
    tmpMask1 = sphdist(pixelRa, pixelDec, ra1, dec1)<cellsize
    tmpMask2 = sphdist(pixelRa, pixelDec, ra2, dec2)<cellsize
    tmpMask4 = sphdist(pixelRa, pixelDec, gra, gdec)<cellsize'''

    tmpMask0 = (abs(ra0-pixelRa)<cellsize)&(abs(dec0-pixelDec)<cellsize)
    tmpMask1 = (abs(ra1-pixelRa)<cellsize)&(abs(dec1-pixelDec)<cellsize)
    tmpMask2 = (abs(ra2-pixelRa)<cellsize)&(abs(dec2-pixelDec)<cellsize)  
    tmpMask4 = (abs(gra-pixelRa)<cellsize)&(abs(gdec-pixelDec)<cellsize)
    print "sum(tmpMask4):", sum(tmpMask4)

    gidRef = gidRef[tmpMask3]
    graRef = graRef[tmpMask3]
    gdecRef = gdecRef[tmpMask3]

    '''minRA = 0.8*min(graRef)
    maxRA = 1.2*max(graRef)
    minDEC = 0.8*min(gdecRef)
    maxDEC = 1.2*max(gdecRef)

    tmpMask0 = (ra0>minRA)&(ra0<maxRA)&(dec0>minDEC)&(dec0<maxDEC)
    tmpMask1 = (ra1>minRA)&(ra1<maxRA)&(dec1>minDEC)&(dec1<maxDEC)
    tmpMask2 =(ra2>minRA)&(ra2<maxRA)&(dec2>minDEC)&(dec2<maxDEC)  
    tmpMask4 = (gra>minRA)&(gra<maxRA)&(gdec>minDEC)&(gdec<maxDEC)'''

    t1= time()
    print "prepare time1", t1 - t0
    obj_id0 = obj_id0[tmpMask0]
    ra0 = ra0[tmpMask0]
    raErr0 = raErr0[tmpMask0]
    dec0 = dec0[tmpMask0]
    decErr0 = decErr0[tmpMask0]
    mjd0 = mjd0[tmpMask0]

    obj_id1 = obj_id1[tmpMask1]
    ra1 = ra1[tmpMask1]
    raErr1 = raErr1[tmpMask1]
    dec1 = dec1[tmpMask1]
    decErr1 = decErr1[tmpMask1]
    mjd1 = mjd1[tmpMask1]

    obj_id2 = obj_id2[tmpMask2]
    ra2 = ra2[tmpMask2]
    raErr2 = raErr2[tmpMask2]
    dec2 = dec2[tmpMask2]
    decErr2 = decErr2[tmpMask2]
    mjd2 = mjd2[tmpMask2]

    gobj_id = gobj_id[tmpMask4]
    gra = gra[tmpMask4]
    graErr = graErr[tmpMask4]
    gdec = gdec[tmpMask4]
    gdecErr = gdecErr[tmpMask4]
    gmjd = gmjd[tmpMask4]
    pixelIndexForObj = pixelIndexForObj[tmpMask4]


    Nref = len(gidRef)

    distTmp = sphdist(pixelRa, pixelDec, graRef, gdecRef)
    angSepMask = np.argsort(distTmp)
    if(Nref>600): Nref = 600
    if(Nref<1): return pixelNo, Nref
    angSepMask = angSepMask[0:Nref]
    objInRadius = gidRef[angSepMask]
    t2= time()
    print "prepare time2:", t2 - t1
    #print "Nref:", Nref, sum(tmpMask0), sum(tmpMask1), sum(tmpMask2), sum(tmpMask3), sum(tmpMask4)
    try:
        cN = 100
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	idx2 = 0
	for idx1 in range(0, Nref):
	    if(idx2>cN-1): break
	    ind0 = obj_id0 == objInRadius[idx1]
	    xtmp0 = mjd0[ind0]
	    xitmp0 = ra0[ind0]
	    etatmp0 = dec0[ind0]

	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1 = mjd1[ind1]
	    xitmp1 = ra1[ind1]
	    etatmp1 = dec1[ind1]	

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2 = mjd2[ind2]
	    xitmp2 = ra2[ind2]
	    etatmp2 = dec2[ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)

	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 
	    if(len(xitmp)<2): continue

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	        raErr2[ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	        decErr2[ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2)):	        
		p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXxi = p[0]*365.0
		slpxi = p[1]

		peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muXeta = peta[0]*365.0
		slpeta = peta[1]

	        gind = gobj_id == objInRadius[idx1]
	        gxtmp = gmjd[gind]
	        gxitmp = gra[gind]
	        getatmp = gdec[gind]

	        offsetXi[idx2] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	        offsetEta[idx2] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
		#print idx2, offsetXi[idx2], offsetEta[idx2]
		idx2 = idx2 + 1

	if(idx2<cN):
	    offsetXi = offsetXi[0:idx2]
	    offsetEta = offsetEta[0:idx2]

    except IndexError:
        return pixelNo, Nref

    t3= time()
    print "prepare time3:", t3 - t2
    indInPixel = pixelIndexForObj == pixelNo
    print "Nref, idx2, cellsize, sum(indInPixel):", Nref, idx1, idx2, cellsize, sum(indInPixel)

    xiStarInPixel  = gra[indInPixel]
    etaStarInPixel = gdec[indInPixel]

    raErrStarInPixel  = graErr[indInPixel]
    decErrStarInPixel = gdecErr[indInPixel]
    starIDInPixel  = gobj_id[indInPixel]
    mjdStarInPixel = gmjd[indInPixel]

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo, Nref
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel


######This function removing the space and magnitude-dependent offset#########
def pixelTasksCombinedDataGAIA5(parameterList):
    # Do not use X-validation, to improve the speed
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, gdra, gddec, gdraREF, gddecREF, pixelIndexForObj, \
	CRA, CDEC = parameterList


    t0= time()

    cellsize = 0.3 #deg

    #tmpMask3 = sphdist(pixelRa, pixelDec, graRef, gdecRef)<cellsize
    tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize) 
    Ntmp = sum(tmpMask3)
    if((Ntmp<110)&(Ntmp>50)): 
	cellsize = 1.0
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)
	Ntmp = sum(tmpMask3)
    if((Ntmp<50)): 
	cellsize = 1.5
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)


    tmpMask0 = (abs(ra0-pixelRa)<cellsize)&(abs(dec0-pixelDec)<cellsize)
    tmpMask1 = (abs(ra1-pixelRa)<cellsize)&(abs(dec1-pixelDec)<cellsize)
    tmpMask2 = (abs(ra2-pixelRa)<cellsize)&(abs(dec2-pixelDec)<cellsize)  
    tmpMask4 = (abs(gra-pixelRa)<cellsize)&(abs(gdec-pixelDec)<cellsize)
    #print "sum(tmpMask4):", sum(tmpMask4)

    gidRef = gidRef[tmpMask3]
    graRef = graRef[tmpMask3]
    gdecRef = gdecRef[tmpMask3]


    t1= time()
    #print "prepare time1", t1 - t0
    obj_id0 = obj_id0[tmpMask0]
    ra0 = ra0[tmpMask0]
    raErr0 = raErr0[tmpMask0]
    dec0 = dec0[tmpMask0]
    decErr0 = decErr0[tmpMask0]
    mjd0 = mjd0[tmpMask0]

    obj_id1 = obj_id1[tmpMask1]
    ra1 = ra1[tmpMask1]
    raErr1 = raErr1[tmpMask1]
    dec1 = dec1[tmpMask1]
    decErr1 = decErr1[tmpMask1]
    mjd1 = mjd1[tmpMask1]

    obj_id2 = obj_id2[tmpMask2]
    ra2 = ra2[tmpMask2]
    raErr2 = raErr2[tmpMask2]
    dec2 = dec2[tmpMask2]
    decErr2 = decErr2[tmpMask2]
    mjd2 = mjd2[tmpMask2]

    gobj_id = gobj_id[tmpMask4]
    gra = gra[tmpMask4]
    graErr = graErr[tmpMask4]
    gdec = gdec[tmpMask4]
    gdecErr = gdecErr[tmpMask4]
    gmjd = gmjd[tmpMask4]
    gdra = gdra[tmpMask4]
    gddec = gddec[tmpMask4]
    pixelIndexForObj = pixelIndexForObj[tmpMask4]


    Nref = len(gidRef)

    distTmp = sphdist(pixelRa, pixelDec, graRef, gdecRef)
    angSepMask = np.argsort(distTmp)
    if(Nref>600): Nref = 600
    if(Nref<1): return pixelNo, Nref
    angSepMask = angSepMask[0:Nref]
    objInRadius = gidRef[angSepMask]
    t2= time()
    #print "prepare time2:", t2 - t1
    #print "Nref:", Nref, sum(tmpMask0), sum(tmpMask1), sum(tmpMask2), sum(tmpMask3), sum(tmpMask4)
    try:
        cN = 100
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	idx2 = 0
	for idx1 in range(0, Nref):
	    if(idx2>cN-1): break
	    ind0 = obj_id0 == objInRadius[idx1]
	    xtmp0 = mjd0[ind0]
	    xitmp0 = ra0[ind0]
	    etatmp0 = dec0[ind0]

	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1 = mjd1[ind1]
	    xitmp1 = ra1[ind1]
	    etatmp1 = dec1[ind1]	

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2 = mjd2[ind2]
	    xitmp2 = ra2[ind2]
	    etatmp2 = dec2[ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)

	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 
	    if(len(xitmp)<2): continue

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	        raErr2[ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	        decErr2[ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2)):	        
		p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXxi = p[0]*365.0
		slpxi = p[1]

		peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muXeta = peta[0]*365.0
		slpeta = peta[1]

	        gind = gobj_id == objInRadius[idx1]
	        gxtmp = gmjd[gind]
	        gxitmp = gra[gind]
	        getatmp = gdec[gind]

	        offsetXi[idx2] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	        offsetEta[idx2] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
		#print idx2, offsetXi[idx2], offsetEta[idx2]
		idx2 = idx2 + 1

	if(idx2<cN):
	    print "pixelNo, Nref, idx2:", pixelNo, Nref, idx2
	    offsetXi = offsetXi[0:idx2]
	    offsetEta = offsetEta[0:idx2]

    except IndexError:
        return pixelNo, Nref

    t3= time()
    #print "prepare time3:", t3 - t2
    indInPixel = pixelIndexForObj == pixelNo
    #print "Nref, idx2, cellsize, sum(indInPixel):", Nref, idx1, idx2, cellsize, sum(indInPixel)

    xiStarInPixel  = gra[indInPixel]
    etaStarInPixel = gdec[indInPixel]

    raErrStarInPixel  = graErr[indInPixel]
    decErrStarInPixel = gdecErr[indInPixel]
    starIDInPixel  = gobj_id[indInPixel]
    mjdStarInPixel = gmjd[indInPixel]

    gdraStarInPixel  = gdra[indInPixel] - gdraREF # remove the offset caused by magnitude
    gddecStarInPixel = gddec[indInPixel] - gddecREF

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    xiStarInPixel += gdraStarInPixel
    etaStarInPixel += gddecStarInPixel

    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo, Nref
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel

######This function removing the space and magnitude-dependent offset#########
def pixelTasksCombinedDataGAIA5DRT(parameterList):
    # Do not use X-validation, to improve the speed
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixelRa, pixelDec, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gidRef, graRef, gdecRef, gobj_id, gra, graErr, \
	gdec, gdecErr, gmjd, gdra, gddec, gdraREF, gddecREF, gpmra, gpmdec, gpmraErr, \
	gpmdecErr, pixelIndexForObj, CRA, CDEC = parameterList


    t0= time()

    cellsize = 0.3 #deg

    #tmpMask3 = sphdist(pixelRa, pixelDec, graRef, gdecRef)<cellsize
    tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize) 
    Ntmp = sum(tmpMask3)
    if((Ntmp<110)&(Ntmp>50)): 
	cellsize = 1.0
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)
	Ntmp = sum(tmpMask3)
    if((Ntmp<50)): 
	cellsize = 1.5
	tmpMask3 = (abs(graRef-pixelRa)<cellsize)&(abs(gdecRef-pixelDec)<cellsize)


    tmpMask0 = (abs(ra0-pixelRa)<cellsize)&(abs(dec0-pixelDec)<cellsize)
    tmpMask1 = (abs(ra1-pixelRa)<cellsize)&(abs(dec1-pixelDec)<cellsize)
    tmpMask2 = (abs(ra2-pixelRa)<cellsize)&(abs(dec2-pixelDec)<cellsize)  
    tmpMask4 = (abs(gra-pixelRa)<cellsize)&(abs(gdec-pixelDec)<cellsize)
    #print "sum(tmpMask4):", sum(tmpMask4)

    gidRef = gidRef[tmpMask3]
    graRef = graRef[tmpMask3]
    gdecRef = gdecRef[tmpMask3]


    t1= time()
    #print "prepare time1", t1 - t0
    obj_id0 = obj_id0[tmpMask0]
    ra0 = ra0[tmpMask0]
    raErr0 = raErr0[tmpMask0]
    dec0 = dec0[tmpMask0]
    decErr0 = decErr0[tmpMask0]
    mjd0 = mjd0[tmpMask0]

    obj_id1 = obj_id1[tmpMask1]
    ra1 = ra1[tmpMask1]
    raErr1 = raErr1[tmpMask1]
    dec1 = dec1[tmpMask1]
    decErr1 = decErr1[tmpMask1]
    mjd1 = mjd1[tmpMask1]

    obj_id2 = obj_id2[tmpMask2]
    ra2 = ra2[tmpMask2]
    raErr2 = raErr2[tmpMask2]
    dec2 = dec2[tmpMask2]
    decErr2 = decErr2[tmpMask2]
    mjd2 = mjd2[tmpMask2]

    gobj_id = gobj_id[tmpMask4]
    gra = gra[tmpMask4]
    graErr = graErr[tmpMask4]
    gdec = gdec[tmpMask4]
    gdecErr = gdecErr[tmpMask4]
    gmjd = gmjd[tmpMask4]
    gdra = gdra[tmpMask4]
    gddec = gddec[tmpMask4]
    gpmra = gpmra[tmpMask4]
    gpmdec = gpmdec[tmpMask4]
    gpmraErr = gpmraErr[tmpMask4]
    gpmdecErr = gpmdecErr[tmpMask4]
    pixelIndexForObj = pixelIndexForObj[tmpMask4]


    Nref = len(gidRef)

    distTmp = sphdist(pixelRa, pixelDec, graRef, gdecRef)
    angSepMask = np.argsort(distTmp)
    if(Nref>900): Nref = 900
    if(Nref<1): return pixelNo, Nref
    angSepMask = angSepMask[0:Nref]
    objInRadius = gidRef[angSepMask]
    t2= time()
    #print "prepare time2:", t2 - t1
    #print "Nref:", Nref, sum(tmpMask0), sum(tmpMask1), sum(tmpMask2), sum(tmpMask3), sum(tmpMask4)
    try:
        cN = 300
	offsetXi = np.zeros(cN)
	offsetEta = np.zeros(cN)
	idx2 = 0
	for idx1 in range(0, Nref):
	    if(idx2>cN-1): break
	    ind0 = obj_id0 == objInRadius[idx1]
	    xtmp0 = mjd0[ind0]
	    xitmp0 = ra0[ind0]
	    etatmp0 = dec0[ind0]

	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1 = mjd1[ind1]
	    xitmp1 = ra1[ind1]
	    etatmp1 = dec1[ind1]	

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2 = mjd2[ind2]
	    xitmp2 = ra2[ind2]
	    etatmp2 = dec2[ind2]
	    #print idx1, np.sum(ind0),  np.sum(ind1), np.sum(ind2)
	    
	    gind = gobj_id == objInRadius[idx1]
	    gxtmp = gmjd[gind]
	    gxitmp = gra[gind]
	    getatmp = gdec[gind]
	    gpmratmp = gpmra[gind]
	    gpmdectmp = gpmdec[gind]
	    #print gpmratmp, gpmdectmp
	    
	    
	    xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))
	    etatmp = np.concatenate((etatmp0, etatmp1, etatmp2)) 
	    if((len(xitmp)<2) or (gpmratmp[0]==np.nan) or (gpmdectmp[0]==np.nan)): continue

	    ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	    yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	        raErr2[ind2]))*1000*3600 

	    ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	    yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	        decErr2[ind2]))*1000*3600
	
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	    if((len(xtmp)>2) & (abs(gpmratmp[0])>0) & (abs(gpmdectmp[0])>0)):
	    	#print gpmratmp, gpmdectmp
	    	p, cov = curve_fit(func2, xtmp*gpmratmp[0]/365.24, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    	muXxi = gpmratmp[0]
	    	slpxi = p[0]
	    	
	    	peta, coveta = curve_fit(func2, xtmp*gpmdectmp[0]/365.24, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	    	muXeta = gpmdectmp[0]
	    	slpeta = peta[0]
	    	
	    	offsetXi[idx2] = ((gxitmp[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * muXxi/365.0 + slpxi))/(1000*3600.0)
	    	offsetEta[idx2] = ((getatmp[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * muXeta/365.0 + slpeta))/(1000*3600.0)
	    	idx2 = idx2 + 1

	if(idx2<cN):
	    print "pixelNo, Nref, idx2:", pixelNo, Nref, idx2
	    offsetXi = offsetXi[0:idx2]
	    offsetEta = offsetEta[0:idx2]

    except IndexError:
        return pixelNo, Nref

    t3= time()
    #print "prepare time3:", t3 - t2
    indInPixel = pixelIndexForObj == pixelNo
    #print "Nref, idx2, cellsize, sum(indInPixel):", Nref, idx1, idx2, cellsize, sum(indInPixel)

    xiStarInPixel  = gra[indInPixel]
    etaStarInPixel = gdec[indInPixel]

    raErrStarInPixel  = graErr[indInPixel]
    decErrStarInPixel = gdecErr[indInPixel]
    starIDInPixel  = gobj_id[indInPixel]
    mjdStarInPixel = gmjd[indInPixel]
    pmraStarInPixel = gpmra[indInPixel]
    pmdecStarInPixel = gpmdec[indInPixel]
    epmraStarInPixel = gpmraErr[indInPixel]
    epmdecStarInPixel = gpmdecErr[indInPixel]

    gdraStarInPixel  = gdra[indInPixel] - gdraREF # remove the offset caused by magnitude
    gddecStarInPixel = gddec[indInPixel] - gddecREF

    medianRaOffsetEpochWise  = np.median(offsetXi)
    medianDecOffsetEpochWise = np.median(offsetEta)
    #update the ra/dec - replace old values
    xiStarInPixel -= medianRaOffsetEpochWise
    etaStarInPixel-= medianDecOffsetEpochWise
    xiStarInPixel += gdraStarInPixel
    etaStarInPixel += gddecStarInPixel

    if((offsetXi.size>1)):
	raErrTmp = np.sqrt((np.pi/2)/(offsetXi.size-1)) * \
	    (0.741*(np.percentile(offsetXi, 75) - np.percentile(offsetXi, 25)))
	decErrTmp = np.sqrt((np.pi/2)/(offsetEta.size-1)) * \
	    (0.741*(np.percentile(offsetEta, 75) - np.percentile(offsetEta, 25)))

	raErrStarInPixel = np.sqrt(raErrStarInPixel**2 + raErrTmp**2)
	decErrStarInPixel = np.sqrt(decErrStarInPixel**2 + decErrTmp**2)

    finalObjIDstar    = starIDInPixel
    #########################(xi, eta) ==> (RA, DEC)#################################
    # transform the (xi, eta)/radians into (RA, DEC)/degree
    finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(xiStarInPixel), np.radians(etaStarInPixel), CRA, CDEC)
    #################################################################################

    if(len(finalObjIDstar)<1):
	return pixelNo, Nref
    else:
	return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, mjdStarInPixel,pmraStarInPixel,pmdecStarInPixel,epmraStarInPixel,epmdecStarInPixel


def pixelTasksCombinedDataCPM(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd , starIDs, starMjds, raStar, decStar, raErrStar, decErrStar, rMagStar, band, x_fix, y_fix, pixAllStar, mjdBreakAt, gnomonic, ra0, dec0 = parameterList

    distTmp = sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)
    '''angSepMask = (distTmp) <= (60.0/60.0)
    if(sum(angSepMask)>600):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:NGAL]
    elif(sum(angSepMask)<100):
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:100]'''

    angSepMask = np.argsort(distTmp)
    angSepMask = angSepMask[0:NGAL]
    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    decFinalGalInRadius = decFinalGal[angSepMask]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	raFinalGalInRadius, decFinalGalInRadius, status = s2t.ds2tp(raFinalGalInRadius, decFinalGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################

    try:
        flag1 = 1
        r1 = time()
        temp = pd.Index(galIDs)
        inde = temp.get_indexer_for(uniqueGalIDinRadius) # all the detections of gal
        inde = inde[inde > -1]
        inde = np.array(inde, dtype = 'i8')
        raGalInRadius  = raGal[inde]
        decGalInRadius = decGal[inde]
	#########################(RA, DEC) ==> (xi, eta)#################################
	if(gnomonic):
	    #print "raGalInRadius, decGalInRadius:", raGalInRadius, decGalInRadius
	    raGalInRadius, decGalInRadius, status = s2t.ds2tp(raGalInRadius, decGalInRadius, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	    #print "Number points in bad transformation:", np.sum(status)
	#################################################################################
        galIDinRadius  = galIDs[inde]
        galMjdInRadius = galMjd[inde]

    except IndexError:
        return pixelNo, flag1, inde.size
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    #print "galIDinRadius:", galIDinRadius
    
    galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    #print "noOfDetectionsOfGal",noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius  
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
                #print "same object now"
        else:
                #make them numpy arrays
		#print "currentGalID, uniqueGalIDinRadius:", currentGalID, uniqueGalIDinRadius[galCounter]
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
                offsetRa  = currentRaGal - raFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetDec = currentDecGal - decFinalGalInRadius[uniqueGalIDinRadius==currentGalID]
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
		#print "offsetRa, offsetDec", offsetRa, offsetDec
                #update counters and variables
	        galCounter+= 1
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]
    #print "time taken to calculate offsets of galaxies from their final ra/dec", time()-t
    #print "galCounter after loop", galCounter	

	
    # select stars inside the pixel
    resIndexInPixel =  pixAllStar == pixelNo
    #print "pixelNo:", pixelNo
    #print "sum of star number in Pixel:", np.sum(resIndexInPixel)
    #uniqueStarIDinPixel = starID[resIndexInPixel]
    #match all detections of stars with those in the pixel
    #indexInPixel = np.in1d(starIDs, uniqueStarIDinPixel)
    #obtain original ra and dec values for stars within pixel
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    #########################(RA, DEC) ==> (xi, eta)#################################
    if(gnomonic):
	#print "raStarInPixel, decStarInPixel:", raStarInPixel, decStarInPixel
	raStarInPixel, decStarInPixel, status = s2t.ds2tp(raStarInPixel, decStarInPixel, ra0, dec0)# transform the (RA, DEC)/degree into (xi, eta)/degree
	#print "Number points in bad transformation:", np.sum(status)
    #################################################################################
    raErrStarInPixel  = raErrStar[resIndexInPixel]
    decErrStarInPixel = decErrStar[resIndexInPixel]
    rMagStarInPixel = rMagStar[resIndexInPixel]
    bandStarInPixel = band[resIndexInPixel]
    x_fixStarInPixel = x_fix[resIndexInPixel]
    y_fixStarInPixel = y_fix[resIndexInPixel]
    #raStarErrInPixel  = np.zeros(np.sum(resIndexInPixel))
    #decStarErrInPixel = np.zeros(np.sum(resIndexInPixel))
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]
    #print 'sum number of resIndexInPixel', np.sum(resIndexInPixel)

    medianRaOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    medianDecOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    #avgMjd = np.zeros(len(mjdBreakAt)+1)

    if (len(mjdBreakAt)>= 1.0):
        for mjdIdx in range(0, len(mjdBreakAt)+1):
	    if(mjdIdx==0):
	        var = mjdBreakAt[mjdIdx]
		mjdIndexInPixel = mjdStarInPixel < var
		mjdIndex = galMjdInRadius < var
	    elif(mjdIdx==len(mjdBreakAt)):
	        var = mjdBreakAt[mjdIdx-1]
		mjdIndexInPixel = mjdStarInPixel >= var
		mjdIndex = galMjdInRadius >= var
	    else:
	        var = mjdBreakAt[mjdIdx]
		mjdIndex = (galMjdInRadius >= mjdBreakAt[mjdIdx-1]) & (galMjdInRadius < var)
                mjdIndexInPixel = (mjdStarInPixel >= mjdBreakAt[mjdIdx-1]) & (mjdStarInPixel < var)
            if (any(mjdIndex)):
                offsetRaValues  = offsetRaArray[mjdIndex]
                offsetDecValues = offsetDecArray[mjdIndex]
                medianRaOffsetEpochWise[mjdIdx]  = np.median(offsetRaValues)
                medianDecOffsetEpochWise[mjdIdx] = np.median(offsetDecValues)
                #update the ra/dec - replace old values
                raStarInPixel[mjdIndexInPixel] -= medianRaOffsetEpochWise[mjdIdx]
                decStarInPixel[mjdIndexInPixel]-= medianDecOffsetEpochWise[mjdIdx]
		#print "np.median(offsetRaValues):", np.median(offsetRaValues)
		if((offsetRaValues.size>1)):
		    raErrTmp = np.sqrt((np.pi/2)/(offsetRaValues.size-1)) * \
		        (0.741*(np.percentile(offsetRaValues, 75) - np.percentile(offsetRaValues, 25)))
		    decErrTmp = np.sqrt((np.pi/2)/(offsetDecValues.size-1)) * \
		        (0.741*(np.percentile(offsetDecValues, 75) - np.percentile(offsetDecValues, 25)))
		    raErrStarInPixel[mjdIndexInPixel] = np.sqrt(raErrStarInPixel[mjdIndexInPixel]**2 + raErrTmp**2) # 
		    decErrStarInPixel[mjdIndexInPixel] = np.sqrt(decErrStarInPixel[mjdIndexInPixel]**2 + decErrTmp**2) # 
                #avgMjd[var1] = (mjdStarInPixel[mjdIndexInPixel].max() - mjdStarInPixel[mjdIndexInPixel].min() )/2

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel

    #########################(xi, eta) ==> (RA, DEC)#################################
    if(gnomonic): # transform the (xi, eta)/radians into (RA, DEC)/degree
	finalRaArrayStar, finalDecArrayStar = s2t.dtp2s(np.radians(finalRaArrayStar), np.radians(finalDecArrayStar), ra0, dec0)
    #################################################################################

    return finalObjIDstar, finalRaArrayStar, finalDecArrayStar, raErrStarInPixel, decErrStarInPixel, \
	mjdStarInPixel,rMagStarInPixel, bandStarInPixel, x_fixStarInPixel, y_fixStarInPixel



def pixelQSOtask(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration

    pixelNo, pixRa, pixDec, objIDs, ra, dec, pmra, epmra, pmdec, epmdec, mr, pixsize = parameterList
	

    cellsize = pixsize #deg

    distTmp = sphdist(pixRa, pixDec, ra, dec)
    angSepMask = (distTmp) <= cellsize

    #angSepMask = np.argsort(distTmp)
    #angSepMask = angSepMask[0:100]

    objInRadius = objIDs[angSepMask]

    pmraInRadius = pmra[angSepMask]
    pmdecInRadius = pmdec[angSepMask]
    epmraInRadius  = epmra[angSepMask]
    epmdecInRadius = epmdec[angSepMask]
    mrInRadius = mr[angSepMask]
    numInRadius = len(pmraInRadius)
    if(numInRadius<10):
	print "num in the pixel:", numInRadius
	return pixelNo, numInRadius

    if (numInRadius>4):
	#print pmraInRadius
	epmraTmp = np.sqrt((np.pi/2)/(numInRadius-1)) * (0.741*(np.percentile(pmraInRadius, 75) - np.percentile(pmraInRadius, 25))) # \
	epmdecTmp = np.sqrt((np.pi/2)/(numInRadius-1)) * (0.741*(np.percentile(pmdecInRadius, 75) - np.percentile(pmdecInRadius, 25)))
	#print epmraTmp, epmdecTmp
        mrapm  = np.median(pmraInRadius)#/epmraTmp
        mdecpm  = np.median(pmdecInRadius)#/epmdecTmp

	epmraInPixel =  (0.741*(np.percentile(pmraInRadius, 75) - np.percentile(pmraInRadius, 25)))#np.sqrt(epmraInRadius**2 + epmraTmp**2)
	epmdecInPixel = (0.741*(np.percentile(pmdecInRadius, 75) - np.percentile(pmdecInRadius, 25)))#np.sqrt(epmdecInRadius**2 + epmdecTmp**2)

    #print pixelNo, mrapm, mdecpm, epmraInPixel, epmdecInPixel, numInRadius
    return pixelNo, pixRa, pixDec, mrapm, mdecpm, epmraInPixel, epmdecInPixel, numInRadius

##################################################################
def QSOpixels(parameterListForPixel, datname):
    dat  = np.zeros(len(parameterListForPixel), dtype=[('pixelNo', 'i4'), ('ra', 'f8'), ('dec', 'f8'), \
	('mrapm', 'f4'), ('mdecpm', 'f4'), ('epmraInPixel', 'f4'), ('epmdecInPixel', 'f4'), ('numInRadius', 'i4')])
    #start workers
    pool = Pool(processes=30)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(pixelQSOtask, parameterListForPixel[0:], chunksize=50)
    noIterated = 0
    for res in iterator:
        if (len(res)==2):
	    print "there was an error with pandas :( ", res
	else:
	    pixelNo_, ra_, dec_, mrapm_, mdecpm_, epmraInPixel_, epmdecInPixel_, numInPixel_ = res
	    #tmp = [int(pixelNo_), (float(mrapm_), float(mdecpm_), float(epmraInPixel_), float(epmdecInPixel_), int(numInPixel_))]
	    dat[noIterated]['pixelNo'] = int(pixelNo_)
	    dat[noIterated]['ra'] = ra_
	    dat[noIterated]['dec'] = dec_
	    dat[noIterated]['mrapm'] = mrapm_
	    dat[noIterated]['mdecpm'] = mdecpm_
	    dat[noIterated]['epmraInPixel'] = epmraInPixel_
	    dat[noIterated]['epmdecInPixel'] = epmdecInPixel_
	    dat[noIterated]['numInRadius'] = int(numInPixel_)
	    noIterated += 1
    pool.terminate()
    dat = dat[0:noIterated]
    np.save(datname, dat)
    print "total pixel nums:", noIterated


##################################################################
def correctedStars(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
    	obj_id = tables.Int64Col(pos=0)
        ra = tables.Float64Col(pos=1)
        dec = tables.Float64Col(pos=2)
	#raErr = tables.Float64Col(pos=3)
	#decErr = tables.Float64Col(pos=4)
        mjd  = tables.Float64Col(pos=3)
        mr = tables.Float64Col(pos=4)
        nObs = tables.Float64Col(pos=5)

    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=100563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=20)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(pixelTasksCombinedData,packParameterList, chunksize=500)
    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==2):
            print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, mjd_, mr_, nObs_ = res
	    noOfPixelsIterated += 1
	    print "noOfPixelsIterated",noOfPixelsIterated
            #so that you donot encounter empty pixels!
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (np.array(objID_)[i]).astype('i8')
		    star['ra']  = (np.array(ra_)[i]).astype('float')
		    star['dec']  = (np.array(dec_)[i]).astype('float')
		    #star['raErr']  = (raErr_[i]).astype('float')
		    #star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (np.array(mjd_)[i]).astype('float')
		    star['mr']  = (np.array(mr_)[i]).astype('float')
		    star['nObs']  = (np.array(nObs_)[i]).astype('float')
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 10000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    #print "no of rows in table we began with", noOfRows
    #terminate the pool of multiprocessors
    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    print "noOfRowsInTable:", noOfRowsInTable
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

##################################################################
def correctedStarsSDSS(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
	obj_id = tables.Int64Col(pos=0)
	ra  = tables.Float64Col(pos=1)
	dec = tables.Float64Col(pos=2)
	raErr = tables.Float64Col(pos=3)
	decErr = tables.Float64Col(pos=4)
	mjd  = tables.Float64Col(pos=5)

    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=30)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    #print "packParameterList[0][17]=", packParameterList[0][17]
    if(packParameterList[0][17][0]==0):
        iterator = pool.imap_unordered(pixelTasksCombinedDataSDSS2,packParameterList, chunksize=100)
    else:
        iterator = pool.imap_unordered(pixelTasksCombinedDataSDSS,packParameterList, chunksize=100)

    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==3):
	    print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, raErr_, decErr_, mjd_ = res
	    noOfPixelsIterated += 1
	    #print "noOfPixelsIterated",noOfPixelsIterated
            #so that you donot encounter empty pixels!
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (objID_[i]).astype('i8')
		    star['ra']  = (ra_[i]).astype('float')
		    star['dec']  = (dec_[i]).astype('float')
		    star['raErr']  = (raErr_[i]).astype('float')
		    star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (mjd_[i]).astype('float')
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 10000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    #print "no of rows in table we began with", noOfRows
    #terminate the pool of multiprocessors
    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    #print "noOfRowsInTable:", noOfRowsInTable
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()


##################################################################
def correctedStarsGAIADRT(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
	obj_id = tables.Int64Col(pos=0)
	ra  = tables.Float64Col(pos=1)
	dec = tables.Float64Col(pos=2)
	raErr = tables.Float64Col(pos=3)
	decErr = tables.Float64Col(pos=4)
	mjd  = tables.Float64Col(pos=5)
	pmra  = tables.Float32Col(pos=6)
	pmdec = tables.Float32Col(pos=7)
	epmra  = tables.Float32Col(pos=8)
	epmdec = tables.Float32Col(pos=9)

	
    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=30)
    ti = time()

    if(packParameterList[0][21][0]==0):
        iterator = pool.imap_unordered(pixelTasksCombinedDataGAIADRT,packParameterList, chunksize=100)

    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==3):
	    print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, raErr_, decErr_, mjd_, pmra_, pmdec_, epmra_, epmdec_ = res
	    noOfPixelsIterated += 1
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (objID_[i]).astype('i8')
		    star['ra']  = (ra_[i]).astype('float')
		    star['dec']  = (dec_[i]).astype('float')
		    star['raErr']  = (raErr_[i]).astype('float')
		    star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (mjd_[i]).astype('float')
		    star['pmra']  = (pmra_[i]).astype('float')
		    star['pmdec']  = (pmdec_[i]).astype('float')
		    star['epmra']  = (epmra_[i]).astype('float')
		    star['epmdec']  = (epmdec_[i]).astype('float')		    
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 10000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    h5file.root.table1.remove()    
    h5file.close()
    
    
 ##################################################################
def correctedStarsGAIADRT2(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
	obj_id = tables.Int64Col(pos=0)
	ra  = tables.Float64Col(pos=1)
	dec = tables.Float64Col(pos=2)
	raErr = tables.Float64Col(pos=3)
	decErr = tables.Float64Col(pos=4)
	mjd  = tables.Float64Col(pos=5)
	pmra  = tables.Float32Col(pos=6)
	pmdec = tables.Float32Col(pos=7)
	epmra  = tables.Float32Col(pos=8)
	epmdec = tables.Float32Col(pos=9)

	
    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=30)
    ti = time()

    iterator = pool.imap_unordered(pixelTasksCombinedDataGAIA5DRT,packParameterList, chunksize=100)

    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==2):
	    print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, raErr_, decErr_, mjd_, pmra_, pmdec_, epmra_, epmdec_ = res
	    noOfPixelsIterated += 1
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (objID_[i]).astype('i8')
		    star['ra']  = (ra_[i]).astype('float')
		    star['dec']  = (dec_[i]).astype('float')
		    star['raErr']  = (raErr_[i]).astype('float')
		    star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (mjd_[i]).astype('float')
		    star['pmra']  = (pmra_[i]).astype('float')
		    star['pmdec']  = (pmdec_[i]).astype('float')
		    star['epmra']  = (epmra_[i]).astype('float')
		    star['epmdec']  = (epmdec_[i]).astype('float')		    
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 10000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    h5file.root.table1.remove()    
    h5file.close()
##################################################################
def correctedStarsGAIA(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
	obj_id = tables.Int64Col(pos=0)
	ra  = tables.Float64Col(pos=1)
	dec = tables.Float64Col(pos=2)
	raErr = tables.Float64Col(pos=3)
	decErr = tables.Float64Col(pos=4)
	mjd  = tables.Float64Col(pos=5)

    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=20)
    ti = time()
    iterator = pool.imap_unordered(pixelTasksCombinedDataGAIA5, packParameterList, chunksize=1000)

    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==2):
	    print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, raErr_, decErr_, mjd_ = res
	    noOfPixelsIterated += 1
	    #print "noOfPixelsIterated",noOfPixelsIterated
            #so that you donot encounter empty pixels!
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (objID_[i]).astype('i8')
		    star['ra']  = (ra_[i]).astype('float')
		    star['dec']  = (dec_[i]).astype('float')
		    star['raErr']  = (raErr_[i]).astype('float')
		    star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (mjd_[i]).astype('float')
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 20000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    #print "no of rows in table we began with", noOfRows
    #terminate the pool of multiprocessors
    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    #print "noOfRowsInTable:", noOfRowsInTable
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()




##############predict each GAIA star, and get the positional offset between predict and orignal position, parallelly ################
def offsetTask(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	CRA, CDEC = packParameterList

    flag = 0  
    ind0 = obj_id0 == uniqueID
    xtmp0 = mjd0[ind0]
    
    if(len(xtmp0)>2):
	ra0 = ra0[ind0]
	raErr0 = raErr0[ind0]
	dec0 = dec0[ind0]
	decErr0 = decErr0[ind0]
	mr0 = mr0[ind0]

	ind1 = obj_id1 == uniqueID
	xtmp1 = mjd1[ind1]
	ra1 = ra1[ind1]
	raErr1 = raErr1[ind1]
	dec1 = dec1[ind1]
	decErr1 = decErr1[ind1]
	if(len(xtmp1)>0): flag = flag + 10	

	ind2 = obj_id2 == uniqueID
	xtmp2 = mjd2[ind2]
	ra2 = ra2[ind2]
	raErr2 = raErr2[ind2]
	dec2 = dec2[ind2]
	decErr2 = decErr2[ind2]
	if(len(xtmp2)>0): flag = flag + 5

	gind = gobj_id == uniqueID
	gxtmp = gmjd[gind]
	gra = gra[gind]
	graErr = graErr[gind]
	gdec = gdec[gind]
	gdecErr = gdecErr[gind]
	if(len(gxtmp)>0): flag = flag + 20

	xitmp = np.concatenate((ra1, ra2, ra0))#   #
	etatmp = np.concatenate((dec1, dec2, dec0))#  # 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0))*1000*3600 # #
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0))#  

	if((len(xtmp)>2)):
	    pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)

	    offsetXi = -1*((gra[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * pXI[0] + pXI[1]))
	    offsetEta = -1*((gdec[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * pETA[0] + pETA[1]))


        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(gra[0]), np.radians(gdec[0]), CRA, CDEC)
        #################################################################################
        return idx, uniqueID, ra, dec, mr0[0], offsetXi, offsetEta, flag
    else:
	return idx, uniqueID


##############predict each GAIA DR2 star, and get the positional offset between predict and orignal position, parallelly ################
def offsetTaskDRT(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, gpmra, gpmdec, \
	CRA, CDEC = packParameterList

    flag = 0  
    ind0 = obj_id0 == uniqueID
    xtmp0 = mjd0[ind0]
    
    if(len(xtmp0)>2):
	ra0 = ra0[ind0]
	raErr0 = raErr0[ind0]
	dec0 = dec0[ind0]
	decErr0 = decErr0[ind0]
	mr0 = mr0[ind0]

	ind1 = obj_id1 == uniqueID
	xtmp1 = mjd1[ind1]
	ra1 = ra1[ind1]
	raErr1 = raErr1[ind1]
	dec1 = dec1[ind1]
	decErr1 = decErr1[ind1]
	if(len(xtmp1)>0): flag = flag + 10	

	ind2 = obj_id2 == uniqueID
	xtmp2 = mjd2[ind2]
	ra2 = ra2[ind2]
	raErr2 = raErr2[ind2]
	dec2 = dec2[ind2]
	decErr2 = decErr2[ind2]
	if(len(xtmp2)>0): flag = flag + 5

	gind = gobj_id == uniqueID
	gxtmp = gmjd[gind]
	gra = gra[gind]
	graErr = graErr[gind]
	gdec = gdec[gind]
	gdecErr = gdecErr[gind]
	gpmra = gpmra[gind]
	gpmdec = gpmdec[gind]
	if(len(gxtmp)>0): flag = flag + 20

	xitmp = np.concatenate((ra1, ra2, ra0))#   #
	etatmp = np.concatenate((dec1, dec2, dec0))#  # 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0))*1000*3600 # #
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0))#  

	if((len(xtmp)>2)):
	    pXI, covXI = curve_fit(func2, xtmp*gpmra[0]/365.24, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    pETA, covETA = curve_fit(func2, xtmp*gpmdec[0]/365.24, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)

	    offsetXi = -1*((gra[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * gpmra[0]/365.24 + pXI[0]))
	    offsetEta = -1*((gdec[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * gpmdec[0]/365.24 + pETA[0]))


        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(gra[0]), np.radians(gdec[0]), CRA, CDEC)
        #################################################################################
        return idx, uniqueID, ra, dec, mr0[0], offsetXi, offsetEta, flag
    else:
	return idx, uniqueID
##################################################################
def poffsetGAIA(packParameterList, poffsetfile):

    poffsetArray = np.zeros(len(packParameterList),  dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), \
	('mr', 'f4'), ('dra', 'f4'), ('ddec', 'f4'), ('flag', 'u2')])
    #start workers
    pool = Pool(processes=20)
    iterator = pool.imap_unordered(offsetTaskDRT, packParameterList, chunksize=500)

    for res in iterator:
        if (len(res)==2):
	    print "there was an error with :( ", res
	else:
            poffsetArray['objID'][res[0]] = res[1]
            poffsetArray['ra'][res[0]] = res[2]
            poffsetArray['dec'][res[0]] = res[3]
            poffsetArray['mr'][res[0]] = res[4]
            poffsetArray['dra'][res[0]] = res[5]
            poffsetArray['ddec'][res[0]] = res[6]
            poffsetArray['flag'][res[0]] = res[7]

    mask = abs(poffsetArray['mr'])>0
    np.save(poffsetfile, poffsetArray[mask])
    pool.terminate()
    return poffsetArray[mask]

##############predict each SDSS star, and get the positional offset between predict and orignal position, parallelly ################
def offsetTask2(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, \
	CRA, CDEC = packParameterList

    flag = 0  
    ind0 = obj_id0 == uniqueID
    xtmp0 = mjd0[ind0]
    
    if(len(xtmp0)>2):
	ra0 = ra0[ind0]
	raErr0 = raErr0[ind0]
	dec0 = dec0[ind0]
	decErr0 = decErr0[ind0]
	mr0 = mr0[ind0]

	ind1 = obj_id1 == uniqueID
	xtmp1 = mjd1[ind1]
	ra1 = ra1[ind1]
	raErr1 = raErr1[ind1]
	dec1 = dec1[ind1]
	decErr1 = decErr1[ind1]
	if(len(xtmp1)>0): flag = flag + 10	



	xitmp = ra0#   #
	etatmp = dec0#  # 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = raErr0*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = decErr0*1000*3600 # #
	
	xtmp = xtmp0  

	if((len(xtmp)>2)):
	    pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)

	    offsetXi = -1*((ra1[0]-np.min(xitmp))*1000*3600 - (xtmp1[0] * pXI[0] + pXI[1]))
	    offsetEta = -1*((dec1[0]-np.min(etatmp))*1000*3600 - (xtmp1[0] * pETA[0] + pETA[1]))


        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
        return idx, uniqueID, ra, dec, mr0[0], offsetXi, offsetEta, flag
    else:
	return idx, uniqueID

##################################################################
def poffsetSDSS(packParameterList, poffsetfile):

    poffsetArray = np.zeros(len(packParameterList),  dtype = [('objID', 'i8'), ('ra', 'f8'), ('dec', 'f8'), \
	('mr', 'f4'), ('dra', 'f4'), ('ddec', 'f4'), ('flag', 'u2')])
    #start workers
    pool = Pool(processes=20)
    iterator = pool.imap_unordered(offsetTask2, packParameterList, chunksize=500)

    for res in iterator:
        if (len(res)==2):
	    print "there was an error with :( ", res
	else:
            poffsetArray['objID'][res[0]] = res[1]
            poffsetArray['ra'][res[0]] = res[2]
            poffsetArray['dec'][res[0]] = res[3]
            poffsetArray['mr'][res[0]] = res[4]
            poffsetArray['dra'][res[0]] = res[5]
            poffsetArray['ddec'][res[0]] = res[6]
            poffsetArray['flag'][res[0]] = res[7]

    mask = abs(poffsetArray['mr'])>0
    np.save(poffsetfile, poffsetArray[mask])
    pool.terminate()
    return poffsetArray[mask]

##################################################################
def correctedStarsCPM(packParameterList, h5correctedStarsFile):
    # table definition
    class Star(tables.IsDescription):
	# in LSD, obj_id is a 64-bit unsigned integer,
	# but pytables cannot index 64-bit unsigned integers
	# so they are saved as 64-bit signed integers
	obj_id = tables.Int64Col(pos=0)
	ra  = tables.Float64Col(pos=1)
	dec = tables.Float64Col(pos=2)
	raErr = tables.Float64Col(pos=3)
	decErr = tables.Float64Col(pos=4)
	mjd  = tables.Float64Col(pos=5)
	mr = tables.Float64Col(pos=6)
	band = tables.StringCol(pos=7, itemsize=1)
	x_fix = tables.Float64Col(pos=8)
	y_fix = tables.Float64Col(pos=9)

    # open a pytable file
    h5file = tables.open_file(h5correctedStarsFile, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    #start workers
    pool = Pool(processes=50)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(pixelTasksCombinedDataCPM,packParameterList, chunksize=100)
    counterForRows = 0
    noOfPixelsIterated = 0
    dat  = []
    for res in iterator:
        if (len(res)==3):
	    print "there was an error with :( ", res
	else:
	    objID_, ra_, dec_, raErr_, decErr_, mjd_, mr_, band_, x_fix_, y_fix_ = res
	    noOfPixelsIterated += 1
	    #print "noOfPixelsIterated",noOfPixelsIterated
            #so that you donot encounter empty pixels!
	    if (objID_.size != 0):
		for i in range(0, objID_.size):    
		    star['obj_id']  = (objID_[i]).astype('i8')
		    star['ra']  = (ra_[i]).astype('float')
		    star['dec']  = (dec_[i]).astype('float')
		    star['raErr']  = (raErr_[i]).astype('float')
		    star['decErr']  = (decErr_[i]).astype('float')
		    star['mjd']  = (mjd_[i]).astype('float')
		    star['mr']  = (mr_[i]).astype('float')
		    star['band']  = (band_[i]).astype('str')
		    star['x_fix']  = (x_fix_[i]).astype('float')
		    star['y_fix']  = (y_fix_[i]).astype('float')
		    star.append()
		    counterForRows += 1
        	if (counterForRows % 10000 == 0): 
            	    table.flush()
	    if (objID_.size ==0):
	 	print "zero array"

    #print "no of rows in table we began with", noOfRows
    #terminate the pool of multiprocessors
    pool.terminate()
    table.flush()
    noOfRowsInTable = table.nrows
    #print "noOfRowsInTable:", noOfRowsInTable
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='correctedStar', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()

def Vpmeq2UVW(Vr,pmra,pmdec,ra,dec,D,Vrerr,pmraerr,pmdecerr,Derr):
    '''#Vpm2UVW converts radial velocities and proper motions to 3D velocities.
    #Parameters:
    #   Vr: radial velocity from observation
    #   pmra, pmdec: proper motion in equatorial coordinates
    #   ra,dec:equatorial coordinates
    #   D: distance to the Sun in pc.
    #Return:
    #   U, V, W: 3D velocity components with respect to LSR
    #REF:Johnson & Soderblom(1987); Klement PhD thesis'''
    velsun=np.array([0,0,0])
    k=4.74057
    ra0=ra*np.pi/180
    dec0=dec*np.pi/180
    cra=np.cos(ra0)
    sra=np.sin(ra0)
    cdec=np.cos(dec0)
    sdec=np.sin(dec0)
    B11=-0.063491*cra*cdec-0.86554 *sra*cdec-0.496799*sdec
    B12= 0.063491*sra-0.86554 *cra+0
    B13= 0.063491*cra*sdec+0.86554 *sra*sdec-0.496799*cdec
    #print B11
    B21= 0.493076*cra*cdec-0.460007*sra*cdec+0.738424*sdec
    B22=-0.493076*sra-0.460007*cra+0
    B23=-0.493076*cra*sdec+0.460007*sra*sdec+0.738424*cdec
    B31=-0.867666*cra*cdec-0.198076*sra*cdec+0.455984*sdec
    B32= 0.867666*sra-0.198076*cra+0
    B33= 0.867666*cra*sdec +0.198076*sra*sdec +0.455984*cdec
    U = velsun[0] + Vr*B11 + k*pmra*D*B12/1000.0 + k*pmdec*D*B13/1000
    V = velsun[1] + Vr*B21 + k*pmra*D*B22/1000.0 + k*pmdec*D*B23/1000
    W = velsun[2] + Vr*B31 + k*pmra*D*B32/1000.0 + k*pmdec*D*B33/1000
    Uerr=np.sqrt(B11**2*Vrerr**2+(k*B12/1000)**2*(D**2*pmraerr**2+pmra**2*Derr**2)+(k*B13/1000)**2*(Derr**2*pmdec**2+pmdecerr**2*D**2))
    Verr=np.sqrt(B21**2*Vrerr**2+(k*B22/1000)**2*(D**2*pmraerr**2+pmra**2*Derr**2)+(k*B23/1000)**2*(Derr**2*pmdec**2+pmdecerr**2*D**2))
    Werr=np.sqrt(B31**2*Vrerr**2+(k*B32/1000)**2*(D**2*pmraerr**2+pmra**2*Derr**2)+(k*B33/1000)**2*(Derr**2*pmdec**2+pmdecerr**2*D**2))
    return U, V, W, Uerr, Verr, Werr

def sigx2(x):
    mx = np.median(x[0])
    return sum(((x[0]-mx)**2)/(len(x[0])-1))
    #return sum(((x[0]-mx)**2 + 4*(x[0]-mx)**2*x[1]**2)/(len(x[0])-1)) # x[0] is UVW, x[1] is the error
def sigxy2(xy):
    mx = np.median(xy[0])
    my = np.median(xy[1])
    return sum(((xy[0]-mx)*(xy[1]-my))/(len(xy[0])-1))
    #return sum(((xy[0]-mx)*(xy[1]-my)+  xy[2]**2 + xy[3]**2)/(len(xy[0])-1))
    #return sum(((xy[0]-mx)*(xy[1]-my)+  (xy[1] - my)**2*xy[2]**2 + (xy[0] - mx)**2*xy[3]**2)/(len(xy[0])-1)) # xy[2] and xy[3] are the errors

def inclined(xy):
    sigxy2 = xy[0]
    #print len(sigxy2)
    sigx2 = xy[1]
    sigy2 = xy[2]
    return 0.5*np.arctan(2*sigxy2/(sigx2 - sigy2))

def ellipTeff(x, U, V, W, Uerr, Verr, Werr, min_x, max_x, binsize, return_empty=False):
    grid = np.arange(min_x, max_x, binsize)
    xC = grid + binsize/2.
    output = np.zeros(xC.size, dtype=[('xC','f4'),('mU','f4'), ('mV','f4'), ('mW','f4'), ('sigU','f4'), ('sigV','f4'), ('sigW','f4'), \
	('sigUV','f4'), ('sigUW','f4'), ('sigVW','f4'), ('UErr','f4'), ('VErr','f4'), ('WErr','f4'), ('sigUErr','f4'), ('sigVErr','f4'), \
	('sigWErr','f4'), ('sigUVErr','f4'), ('sigUWErr','f4'), ('sigVWErr','f4'), ('lv','f4'), ('alpha','f4'), ('lvErr','f4'), ('alphaErr','f4'), ('cnt','i4')])
    Nsample = 300
    for i in np.arange(xC.size):
        xbin = np.where((x > grid[i]) & (x <= grid[i]+binsize))
        count = len(xbin[0])
        if count > 4:
	    output['xC'][i] = xC[i]
    	    '''output['mU'][i]=np.median(U[xbin])
    	    output['mV'][i]=np.median(V[xbin])
    	    output['mW'][i]=np.median(W[xbin])
	    #print mU[i]
    	    sigU2=sigx2([U[xbin], Uerr[xbin]])
    	    sigV2=sigx2([V[xbin], Verr[xbin]])
    	    sigW2=sigx2([W[xbin], Werr[xbin]])
    	    sigUV2=sigxy2([U[xbin], V[xbin], Uerr[xbin], Verr[xbin]])
    	    sigUW2=sigxy2([U[xbin], W[xbin], Uerr[xbin], Werr[xbin]])
    	    sigVW2=sigxy2([V[xbin], W[xbin], Verr[xbin], Werr[xbin]])
    	    output['sigU'][i]=np.sqrt(sigU2)
    	    output['sigV'][i]=np.sqrt(sigV2)
    	    output['sigW'][i]=np.sqrt(sigW2)
    	    output['sigUV'][i]=np.sign(sigUV2)*np.sqrt(abs(sigUV2))
    	    output['sigUW'][i]=np.sign(sigUW2)*np.sqrt(abs(sigUW2))
    	    output['sigVW'][i]=np.sign(sigVW2)*np.sqrt(abs(sigVW2))'''

	    tmpidx = bt.subsample_indexes(U[xbin], n_samples=Nsample, size=0.7)	    
	    tmpU  = np.array([np.median(U[xbin][indexes]) for indexes in tmpidx])
	    output['mU'][i]= np.median(tmpU)
	    output['UErr'][i] = (0.741*(np.percentile(tmpU, 75) - np.percentile(tmpU, 25)))

	    #tmpidx = bt.subsample_indexes(V[xbin], n_samples=Nsample, size=0.6)	    
	    tmpV  = np.array([np.median(V[xbin][indexes]) for indexes in tmpidx])
	    output['mV'][i]= np.median(tmpV)
	    output['VErr'][i] = (0.741*(np.percentile(tmpV, 75) - np.percentile(tmpV, 25)))

	    #tmpidx = bt.subsample_indexes(W[xbin], n_samples=Nsample, size=0.6)	    
	    tmpW  = np.array([np.median(W[xbin][indexes]) for indexes in tmpidx])
	    output['mW'][i]= np.median(tmpW)
	    output['WErr'][i] = (0.741*(np.percentile(tmpW, 75) - np.percentile(tmpW, 25)))

	    #tmpidx = bt.subsample_indexes(np.array([U[xbin], Uerr[xbin]]), n_samples=Nsample, size=0.6)
	    tmpSigU2 = np.array([U[xbin], Uerr[xbin]])
	    #print  "U[xbin][0:10], Uerr[xbin][0:10]:", U[xbin][0:10], Uerr[xbin][0:10]
	    #print "len(tmpSigU2):", len(tmpSigU2.T), tmpSigU2[0][0:10]	    
	    tmpSigU2  = np.array([sigx2((tmpSigU2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigU = np.sqrt(tmpSigU2)
	    output['sigU'][i]= np.median(tmpSigU)
	    output['sigUErr'][i] = (0.741*(np.percentile(tmpSigU, 75) - np.percentile(tmpSigU, 25)))
	    #print output['sigU'][i], output['sigUErr'][i], tmpSigU[0:30]

	    tmpSigV2 = np.array([V[xbin], Verr[xbin]]) 	    
	    tmpSigV2  = np.array([sigx2((tmpSigV2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigV = np.sqrt(tmpSigV2)
	    output['sigV'][i]= np.median(tmpSigV)
	    output['sigVErr'][i] = (0.741*(np.percentile(tmpSigV, 75) - np.percentile(tmpSigV, 25)))

	    tmpSigW2 = np.array([W[xbin], Werr[xbin]]) 	    
	    tmpSigW2  = np.array([sigx2((tmpSigW2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigW = np.sqrt(tmpSigW2)
	    output['sigW'][i]= np.median(tmpSigW)
	    output['sigWErr'][i] = (0.741*(np.percentile(tmpSigW, 75) - np.percentile(tmpSigW, 25)))

	    tmpSigUV2 = np.array([U[xbin], V[xbin], Uerr[xbin], Verr[xbin]])	    
	    tmpSigUV2  = np.array([sigxy2((tmpSigUV2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigUV = np.sign(tmpSigUV2)*np.sqrt(abs(tmpSigUV2))
	    output['sigUV'][i]= np.median(tmpSigUV)
	    output['sigUVErr'][i] = (0.741*(np.percentile(tmpSigUV, 75) - np.percentile(tmpSigUV, 25)))

	    tmpSigUW2 = np.array([U[xbin], W[xbin], Uerr[xbin], Werr[xbin]])	    
	    tmpSigUW2  = np.array([sigxy2((tmpSigUW2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigUW = np.sign(tmpSigUW2)*np.sqrt(abs(tmpSigUW2))
	    output['sigUW'][i]= np.median(tmpSigUW)
	    output['sigUWErr'][i] = (0.741*(np.percentile(tmpSigUW, 75) - np.percentile(tmpSigUW, 25)))

	    tmpSigVW2 = np.array([V[xbin], W[xbin], Verr[xbin], Werr[xbin]])	    
	    tmpSigVW2  = np.array([sigxy2((tmpSigVW2.T[indexes]).T) for indexes in tmpidx])
	    tmpSigVW = np.sign(tmpSigVW2)*np.sqrt(abs(tmpSigVW2))
	    output['sigVW'][i]= np.median(tmpSigVW)
	    output['sigVWErr'][i] = (0.741*(np.percentile(tmpSigVW, 75) - np.percentile(tmpSigVW, 25)))

	    tmpidx = bt.subsample_indexes(tmpSigUV2, n_samples=100, size=0.7)	
	    tmplv = np.array([tmpSigUV2, tmpSigU2, tmpSigV2])	 
	    print len(tmplv.T), tmplv[0]
	    tmplv  = np.array([inclined((tmplv.T[indexes]).T) for indexes in tmpidx])*180.0/np.pi
	    output['lv'][i]= np.median(tmplv)
	    output['lvErr'][i] = (0.741*(np.percentile(tmplv, 75) - np.percentile(tmplv, 25)))

	    tmpAlpha = np.array([tmpSigUW2, tmpSigU2, tmpSigW2])	    
	    tmpAlpha  = np.array([inclined((tmpAlpha.T[indexes]).T) for indexes in tmpidx])*180.0/np.pi
	    output['alpha'][i]= np.median(tmpAlpha)
	    output['alphaErr'][i] = (0.741*(np.percentile(tmpAlpha, 75) - np.percentile(tmpAlpha, 25)))

	    '''output['sigUErr'][i] = np.mean(0.5*(bt.ci((U[xbin], Uerr[xbin]), sigx2, 0.157*2, \
		300, 'bca', 'errorbar', 0.001))/np.sqrt(sigU2))
	    output['sigVErr'][i] = np.mean(0.5*(bt.ci((V[xbin], Verr[xbin]), sigx2, 0.157*2, \
		300, 'bca', 'errorbar', 0.001))/np.sqrt(sigV2))
	    output['sigWErr'][i] = np.mean(0.5*(bt.ci((W[xbin], Werr[xbin]), sigx2, 0.157*2, \
		300, 'bca', 'errorbar', 0.001))/np.sqrt(sigW2))
	    output['sigUVErr'][i] = np.mean(0.5*(bt.ci(((U[xbin], V[xbin], Uerr[xbin], Verr[xbin])), \
		sigxy2, 0.157*2, 300, 'bca', 'errorbar', 0.001))/np.sqrt(abs(sigUV2)))
	    output['sigUWErr'][i] = np.mean(0.5*(bt.ci((U[xbin], W[xbin], Uerr[xbin], Werr[xbin]), \
		sigxy2, 0.157*2, 300, 'bca', 'errorbar', 0.001))/np.sqrt(abs(sigUW2)))
	    output['sigVWErr'][i] = np.mean(0.5*(bt.ci((V[xbin], W[xbin], Verr[xbin], Werr[xbin]), \
		sigxy2, 0.157*2, 300, 'bca', 'errorbar', 0.001))/np.sqrt(abs(sigVW2)))'''
        output['cnt'][i] = count
    if return_empty:
        return output
    else:
        index = np.where( output['cnt'] > 0 )
        output = output[index]
        return output
