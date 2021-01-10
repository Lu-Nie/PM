#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 00:00:02 2020

@author: njl
"""

import numpy as np
import tables
import matplotlib as mpl
import pylab as pl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from esutil.coords import sphdist
from astropy.time import Time
from time import time
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import ICRS, Galactic,SkyCoord,Distance,FK5
from astropy import units as u
import sys
import spherical_to_tangential as s2t
from scipy.optimize import curve_fit
import numpy.ma as ma
import esutil.stat
import useful
import matplotlib.cm as cm
import healpy as hp
import random
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
mpl.rc('font', family='serif', serif='cm10')
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import emcee
import corner

def func(x, a, b):
    return a*x + b

def func2(x, b): # a is included in x
    return x + b
#************************updata by Nie*************************
xtmp = []
rMagtmp = []
R2XI = []
chi2XI = []
muXI = []
muErrXI = []
muXxi = []
muErrXxi = []
R2ETA = []
chi2ETA = []
muETA = []
muErrETA = []
muXeta = []
muErrXeta = []
xpm = []

idx=index
uniqueID=uniqueID
obj_id0=obj_id0
ra0=ra0
dec0=dec0
mjd0=mjd0
mr0=mr0
obj_id1=obj_id1
ra1=ra1
raErr1=raErr1
dec1=dec1
decErr1=decErr1
mjd1=mjd1
obj_id2=obj_id2
ra2=ra2
raErr2=raErr2
dec2=dec2
decErr2=decErr2
mjd2=mjd2
pflag=True
chunkNo=0
CRA = CRA
CDEC = CDE
root = root

flag = 0  #record the PS1 
ind0 = obj_id0 == uniqueID
    #print 'ind0', ind0
if(np.sum(ind0)>0):
	#xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
    xtmp0 = mjd0[ind0]
	#xitmp0 = ra0[ind0]
	#etatmp0 = dec0[ind0]
	#CRA = np.median(ra0[ind0])  # Here is important, CRA, CDEC should keep the same as previous
 	#CDEC = np.median(dec0[ind0])
xitmp0, etatmp0, status = s2t.ds2tp(ra0[ind0], dec0[ind0], CRA, CDEC)
print (xitmp0-min(xitmp0))*3600000, (etatmp0-min(etatmp0))*3600000
ind1 = obj_id1 == uniqueID
	#xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
xtmp1 = mjd1[ind1]
xitmp1, etatmp1, status = s2t.ds2tp(ra1[ind1], dec1[ind1], CRA, CDEC)
if(np.sum(ind1)>0): flag = flag + 10

        #print ra0[ind0], dec0[ind0], ra1[ind0], dec1[ind0]	

ind2 = obj_id2 == uniqueID
	#xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
xtmp2 = mjd2[ind2]
xitmp2, etatmp2, status = s2t.ds2tp(ra2[ind2], dec2[ind2], CRA, CDEC)
if(np.sum(ind2)>0): flag = flag + 5

xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))# 
etatmp = np.concatenate((etatmp0, etatmp1, etatmp2))# 

ytmpXI = (xitmp-np.min(xitmp))*1000*3600
#yErrtmpXI = 0.000000001+np.concatenate((raErr1[ind1], \
#	   raErr2[ind2]))*1000*3600 

ytmpETA = (etatmp-np.min(etatmp))*1000*3600
#yErrtmpETA = 0.000000001+np.concatenate((decErr1[ind1], \
#	   decErr2[ind2]))*1000*3600

global xtmp
xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 
    #global xtmp
	#print "yErrtmpXI", yErrtmpXI, yErrtmpETA, xtmp
	#if(pflag):
global rMagtmp
global R2XI
global chi2XI
global muXI
global muErrXI
global muXxi
global muErrXxi
global R2ETA
global chi2ETA
global muETA
global muErrETA
global muXeta
global muErrXeta
global xpm
yErrtmpXI = np.ones([2])
yErrtmpETA = np.ones([2])
print 'len', len(xtmp)
if((len(xtmp)>0)):
	rMagtmp = mr0[ind0][0]
	print(1)
	num0 = np.sum(ind0)	
		#print "starIdx, num0, chunkNo:", idx, num0,chunkNo
	xpm = np.zeros([6, num0])
	
	
	for idxX in range(0, num0): # X means cross-validation
		#yErrtmpXI[idxX] = yErrtmpXI[idxX] + 9999.0
		#yErrtmpETA[idxX] = yErrtmpETA[idxX] + 9999.0
			#print yErrtmpXI
		pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365
		xpm[2, idxX] = (ytmpXI[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		xpm[3, idxX] = pXeta[0]*365
		xpm[4, idxX] = np.sqrt(np.diag(covXeta))[0]*365
		xpm[5, idxX] = (ytmpETA[idxX] - (pXeta[0]*xtmp[idxX] + pXeta[1]))
		#yErrtmpXI[idxX] = np.array(yErrtmpXI)[idxX] - 9999.0
        yErrtmpXI[idxX] = np.array(yErrtmpXI)[idxX]
        #yErrtmpETA[idxX] = yErrtmpETA[idxX] - 9999.0
        yErrtmpETA[idxX] = np.array(yErrtmpETA)[idxX]
		#print xtmp, ytmpXI, yErrtmpXI
	muXxi = np.median(xpm[0])
	muErrXxi = (0.741*(np.percentile(xpm[0], 75) - np.percentile(xpm[0], 25)))
	muXeta = np.median(xpm[3])
	muErrXeta = (0.741*(np.percentile(xpm[3], 75) - np.percentile(xpm[3], 25)))

	print xtmp, ytmpXI, yErrtmpXI
	pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	muXI = pXI[0]*365
	print muXI
	muErrXI = np.sqrt(np.diag(covXI))[0]*365
	variance = np.var(ytmpXI)
	residuals = np.var(ytmpXI - pXI[0]*xtmp + pXI[1])
	Rtmp = residuals/variance		
	R2XI = np.round(np.abs(1-Rtmp), decimals=2) 
	chi2XI =np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-1) # we should use the reduced chi2
		#print np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/5.0

	pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	muETA = pETA[0]*365
	muErrETA = np.sqrt(np.diag(covETA))[0]*365
	variance = np.var(ytmpETA)
	residuals = np.var(ytmpETA - pETA[0]*xtmp + pETA[1])
	Rtmp = residuals/variance		
	R2ETA = np.round(np.abs(1-Rtmp), decimals=2) 
	chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-1)# we should use the reduced chi2
	if(True):
		print('1')#pXI
		p = pETA#pXI
		mu = muETA#muXI
		muErr = muErrETA#muErrXI
		yErrtmp = yErrtmpETA#yErrtmpXI
		ytmp = ytmpETA#ytmpXI
		muX = muXeta#muXxi
		muErrX = muErrXeta#muErrXxi
		R2 = R2ETA#R2XI
		chi2 = chi2ETA#chi2XI

		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		#ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		ax.errorbar(xtmp, ytmp, [0,0], fmt='.', color='red')
			#ax.plot(xtmp, xtmp*(muX/365)+(np.median(ytmp)-muX/365*np.median(xtmp)), 'b--')
		ax.text(np.min(xtmp)+50, np.max(ytmp)*0.72,'PM = %0.2f'% ytmp[0], fontsize=10)
		ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,r'$\chi^2 = %0.2f$'% chi2, fontsize=10)
		ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,r'$R^2 = %0.2f$'% R2, fontsize=10)
		ax.plot(xtmp, xtmp*p[0]+p[1], 'b-')
		ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
		ax.set_ylabel(r'$\delta$(mas)')
		ax.set_xlabel('MJD')
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, r"$objID = %d$" %uniqueID, fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, r"$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
		
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "mu = %0.2f" % mu +r"+-%0.2f" \
		    	    %muErr + ", flag = %d" %(flag), fontsize=10)

		ax.set_title('Proper Motions with M31', fontsize=12)
		figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
		plt.savefig(root+"figs/%s_Mock.png" %figname)
			#plt.show()
print uniqueID, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, \
	    chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag, ra0[ind0][0], dec0[ind0][0]
    #else:
	#return idx, uniqueID
'''ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)'''