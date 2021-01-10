#!/usr/bin/env python
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
#***************************************************************
def func(x, a, b):
    return a*x + b

def func2(x, b): # a is included in x
    return x + b
    
##############histogram medianDeltaRA for different patches in one MJD################
def histMedianResiMJD(medianMJD, medianRAresidual, stdRAresidual, medianRADEC, noGalaxyMJD, root):
    ind = (np.sum(np.abs(medianRAresidual)>0, 0)>40)
    mjdN = np.sum(ind)
    medianMJD = medianMJD[:,ind]
    medianRAresidual =  medianRAresidual[:,ind]
    #stdRAresidual = stdRAresidual[:,ind]
    fig, axs = plt.subplots(nrows=mjdN, ncols=1, sharex=True, figsize=(4.0,8.0))
    plt.subplots_adjust(left=0.19, bottom=0.09, right=0.99, top=0.99, wspace=0, hspace=0)
    for mjdIdx in range(0, mjdN):	
        ax = axs[mjdIdx]
        tmp = medianRAresidual[:,mjdIdx]
        cri = (tmp>-200)&(tmp<200)#&(np.abs(tmp)>0)
        tmp = tmp[cri]
        Ntmp = np.sum(cri)
        print np.median(tmp)
        if(Ntmp>100): Nbin = 40
        if(Ntmp<100): Nbin = 20
        histo = ax.hist(tmp, int(Nbin), [-20, 20],color='k')
        yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
        #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
        ax.plot([0, 0], [0, 9000], '--', markersize=2.5, color='y')
        ax.set_ylim([0, np.max(histo[0])*1.1])	
        ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
        ax.set_ylabel('$N$')
        ax.annotate('MJD:'+ str(round(np.median(medianMJD[cri,mjdIdx]), 2)), xy=(-15.0, yscale*0.8), \
	    xycoords='data',xytext=(10, 0),color='k', textcoords='offset points', \
            size=8, va="center")

    plt.xlabel('$\Delta RA$/mas')
    #plt.ylabel('N')
    plt.xlim(-20,20) 
    #plt.xticks(np.array([-80, -40, 0, 40, 80]), fontsize='small')
    #plt.show()
    plt.savefig(root+"residualRA_MJD_hist_All.png") 


##############histogram medianDeltaRA for different patches in one rMAG################
def histMedianResiMAG(medianMag, medianRAresiMag, stdRAresiMag, medianRADEC, noGalaxyMag, root):
    fig, axs = plt.subplots(nrows=len(mrBreakAt)+1, ncols=1, sharex=True, figsize=(4.0,8.0))
    plt.subplots_adjust(left=0.19, bottom=0.09, right=0.99, top=0.99, wspace=0., hspace=0.)
    for mrIdx in range(0, len(mrBreakAt)+1):	
        ax = axs[mrIdx]
        tmp = medianRAresiMag[:,mrIdx]
        cri = (tmp>-200)&(tmp<200)&(np.abs(tmp)>0)
        tmp = tmp[cri]
        Ntmp = np.sum(cri)
        print np.median(tmp)
        if(Ntmp>100): Nbin = 50
        if(Ntmp<100): Nbin = 20
        histo = ax.hist(tmp, int(Nbin), [-5, 5],color='k')
        yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
        #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
        ax.plot([0, 0], [0, 9000], '--', markersize=2.5, color='y')
        ax.set_ylim([0, np.max(histo[0])*1.1])	
        ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
        ax.set_ylabel('$N$')
        ax.annotate('rMAG:'+ str(round(np.median(medianMag[cri,mrIdx]), 2)), xy=(-4.0, yscale*0.8), \
	    xycoords='data',xytext=(10, 0),color='k', textcoords='offset points', \
            size=8, va="center")

    plt.xlabel('$\Delta RA$/mas')
    #plt.ylabel('N')
    plt.xlim(-5,5) 
    #plt.xticks(np.array([-80, -40, 0, 40, 80]), fontsize='small')
    plt.savefig(root+"residualRA_Mag_hist_All.png") 

##############plot DetaRA v.s. epoch(MJD) for different patches in one panel################
def scatterMedianResiMJD1(medianMJD, medianRAresidual, stdRAresidual, medianRADEC, noGalaxyMJD, root):
    ind = (np.sum(np.abs(medianRAresidual)>0, 0)>40)
    mjdN = np.sum(ind)
    medianMJD = medianMJD[:,ind]
    medianRAresidual =  medianRAresidual[:,ind]
    noGalaxyMJD = noGalaxyMJD[:,ind]
    fig, axs = plt.subplots(nrows=mjdN, ncols=2, sharex=True, figsize=(6.0,8.0))
    plt.subplots_adjust(left=0.1, bottom=0.09, right=0.99, top=0.94, wspace=0., hspace=0.0)
    print axs.ravel()[0]
    for mjdIdx in range(0, mjdN):	
        tmp = medianRAresidual[:,mjdIdx]
        cri = (tmp>-200)&(tmp<200)&(np.abs(tmp)>0)&(medianRADEC[:,0]<300)&(medianRADEC[:,0]>30)&(noGalaxyMJD[:,mjdIdx]>100)
        tmp = tmp[cri]
        resiTMP = medianRAresidual[cri,mjdIdx]
        ax = axs.ravel()[mjdIdx*2]
        img = ax.scatter(medianRADEC[cri,0], medianRADEC[cri,1], c=resiTMP, cmap=plt.cm.jet, marker='o')
        if(mjdIdx==0):
    	    xscale = (np.max(medianRADEC[cri,0]) - np.min(medianRADEC[cri,0]))
    	    yscale = (np.max(medianRADEC[cri,1]) - np.min(medianRADEC[cri,1]))
	    xlim = np.array([np.min(medianRADEC[cri,0]), np.max(medianRADEC[cri,0])])
	    ylim = np.array([np.min(medianRADEC[cri,1]), np.max(medianRADEC[cri,1])]) 
	    cscale = (np.max(resiTMP) - np.min(resiTMP))
	    clim = np.array([np.min(resiTMP), np.max(resiTMP)])
        #ax.annotate('MJD:'+ str(round(np.median(medianMJD[cri,mjdIdx]), 2)), xy=(20, 28), \
        #	xycoords='data',xytext=(10, 0),color='k', textcoords='offset points', \
        #	size=8, va="center")
        ax.set_xlim([xlim[0]-xscale/20, xlim[1]+xscale/20])
        ax.set_ylim([ylim[0]-yscale/10, ylim[1]+yscale/10])
        ax.set_ylabel('DEC')
        ax.set_yticks([int(ylim[0]+yscale/10), int(ylim[0]+yscale/2), int(ylim[1]-yscale/10)])    
        div = make_axes_locatable(ax)
        cax = div.append_axes("top", size="5%", pad=-0.054)
        print np.max(resiTMP), np.min(resiTMP)
        cbar = plt.colorbar(img, cax=cax, orientation='horizontal')
        cax.set_visible(False)
        cbar.set_clim([-10,10])
        #cbar.set_ticks([-7,0,7])
        if((mjdIdx*2)>(mjdN*2-3)): 
	    ax.set_xlabel('RA')
	    #ax.set_xticks([0, 20, 40, 60, 80])
	    ax.set_xticks([xlim[0]+xscale/10, xlim[0]+xscale/2, xlim[1]-xscale/10])
        if(mjdIdx==0):
	    cax.xaxis.set_ticks_position("top")
	    cax.set_visible(True)
    
        ax = axs.ravel()[mjdIdx*2+1]
        stdTMP = stdRAresidual[cri,mjdIdx]
        img1 = ax.scatter(medianRADEC[cri,0], medianRADEC[cri,1], c=stdTMP, marker='o')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_xlim([xlim[0]-xscale/20, xlim[1]+xscale/20])
        ax.set_ylim([ylim[0]-yscale/10, ylim[1]+yscale/10])
        #print np.min(medianRADEC[cri,1]), np.max(medianRADEC[cri,1]), xscale
        div = make_axes_locatable(ax)
        cax = div.append_axes("top", size="5%", pad=-0.054)
        print np.max(stdTMP), np.min(stdTMP)
        cscale = (np.max(stdTMP) - np.min(stdTMP))
        cbar = plt.colorbar(img1, cax=cax, orientation='horizontal')
        cax.set_visible(False)
        cbar.set_clim(vmin=0, vmax=5.2)
        #cbar.set_ticks([0,1.0,2.0, 3.0, 4.0, 5.0])
        if((mjdIdx*2 + 1)>mjdN*2-3): 
	    ax.set_xlabel('RA')
	    ax.set_xticks([xlim[0]+xscale/10, xlim[0]+xscale/2, xlim[1]-xscale/10])
        if(mjdIdx==0):
	    cax.set_visible(True)
	    cax.xaxis.set_ticks_position("top")
       
    #plt.ylim(-8.0,7.)
    #plt.xticks(np.array([0, 40, 80]), fontsize='small')
    #plt.yticks(np.array([-6, -3, 0, 3, 6]), fontsize='small')
    #plt.title('Galaxies distribution') 
    #plt.grid(True)
    #plt.colorbar()
    plt.savefig(root+"distri_galaxies.png")

##############plot DetaRA v.s. epoch(MJD) for different patches in one panel################
def scatterMedianResiMJD(medianMJD, medianRAresidual, stdRAresidual, medianRADEC, noGalaxyMJD, root):
    Years = ["2010", "2011", "2012", "2013", "2014", "2015", "2016"]
    xmin, xmax, ymin, ymax = 48, 131, -1, 24
    ind = (np.sum(np.abs(medianRAresidual)>0, 0)>50)
    print medianMJD[medianMJD.T[0]>0], np.sum(ind)
    mjdN = np.sum(ind) - 1# Here should not -1
    medianMJD = medianMJD[:,ind]
    medianRAresidual =  medianRAresidual[:,ind]
    stdRAresidual = stdRAresidual[:,ind]*np.sqrt(noGalaxyMJD[:,ind])
    noGalaxyMJD = noGalaxyMJD[:,ind]

    #############Galactic Disk###############
    lDisk = np.arange(0, 360, 1.0)
    bDisk0 = np.zeros(len(lDisk)) -20
    bDisk1 = np.zeros(len(lDisk)) + 0
    bDisk2 = np.zeros(len(lDisk)) + 20
    radec0 = SkyCoord("galactic", l=lDisk, b=bDisk0, unit=(u.degree, u.degree)).fk5
    radec1 = SkyCoord("galactic", l=lDisk, b=bDisk1, unit=(u.degree, u.degree)).fk5
    radec2 = SkyCoord("galactic", l=lDisk, b=bDisk2, unit=(u.degree, u.degree)).fk5
    #print "glactic disk", np.max(radec1.ra.value), np.min(radec1.ra.value), np.max(radec1.dec.value), np.min(radec1.dec.value)
    #print "glactic disk", radec1.ra.value[(radec1.dec.value<25)&(radec1.dec.value>10)], radec1.dec.value[(radec1.dec.value<25)&(radec1.dec.value>10)],
    fig, axs = plt.subplots(nrows=mjdN, ncols=3, sharex=True, figsize=(9.0,2.*(mjdN)))
    #plt.subplots_adjust(left=0.1, bottom=0.09, right=0.99, top=0.9, wspace=0., hspace=0.0)
    plt.subplots_adjust(left=0.1, bottom=0.06, right=0.99, top=0.92, wspace=0., hspace=0.0)
    #print axs.ravel()[0]
    for mjdIdx in range(0, mjdN):	
        tmp = medianRAresidual[:,mjdIdx]
	print np.sum(noGalaxyMJD[:,mjdIdx]<50)
        cri = (tmp>-500)&(tmp<500)&(np.abs(tmp)>0)&(medianRADEC[:,0]<300)&(medianRADEC[:,0]>30)&(noGalaxyMJD[:,mjdIdx]>10)# & (medianMJD.T[0]>51500) & (medianMJD.T[0]<51600)
        #tmp = tmp[cri]
 
        if(mjdIdx==0):
    	    xscale = (np.max(medianRADEC[cri,0]) - np.min(medianRADEC[cri,0]))
    	    yscale = (np.max(medianRADEC[cri,1]) - np.min(medianRADEC[cri,1]))
	    xlim = np.array([(np.min(medianRADEC[cri,0])), (np.max(medianRADEC[cri,0]))])
	    ylim = np.array([(np.min(medianRADEC[cri,1])), (np.max(medianRADEC[cri,1]))]) 
	    #xlim = np.array([248.4, 252.5])
	    #ylim = np.array([34.4, 38.5]) 
	    #cscale0 = int(np.max(resiTMP) - np.min(resiTMP))
	    #clim0 = np.array([int(np.min(resiTMP)), int(np.max(resiTMP))])
	    cax0 = fig.add_axes([0.1, 0.925, 0.2945, 0.016])
	    cmap0 = mpl.cm.jet  #cool
	    norm0 = mpl.colors.Normalize(vmin=-10, vmax=10)
	    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
	    cax0.xaxis.set_ticks_position("top")
	    cbar0.set_ticks([-8,-4, 0,4, 8])
	    cbar0.set_label(r'\textbf{Median Offsets in $\delta$ (mas)}', fontsize='small',fontname='Arial')
	    cbar0.ax.xaxis.set_label_position('top')
	    cax1 = fig.add_axes([0.3963, 0.925, 0.2955, 0.016])
	    cmap1 = mpl.cm.jet  #cool
	    norm1 = mpl.colors.Normalize(vmin=0, vmax=120)#(vmin=0, vmax=7.5)
	    cbar1 = mpl.colorbar.ColorbarBase(cax1, cmap=cmap1,norm=norm1,orientation='horizontal')
	    cax1.xaxis.set_ticks_position("top")
	    cbar1.set_ticks([20,50,80,100])#([1,3,5,7])
	    #cbar1.set_label(r'\textbf{RMS of the Offsets(DEC/mas)}', fontsize='small',fontname='Arial')
	    cbar1.set_label(r'\textbf{RMS of the Offsets in $\delta$ (mas)}', fontsize='small',fontname='Arial')
	    cbar1.ax.xaxis.set_label_position('top')
	    cax2 = fig.add_axes([0.6948, 0.925, 0.295, 0.016])
	    cmap2 = mpl.cm.jet  #cool
	    norm2 = mpl.colors.Normalize(vmin=50, vmax=1500)#(vmin=20, vmax=180)#
	    cbar2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmap2,norm=norm2,orientation='horizontal')
	    cax2.xaxis.set_ticks_position("top")
	    cbar2.set_ticks([300,600,900,1200])#([50,100,150]) #
	    cbar2.set_label(r"\textbf{Galaxies Number in Patches}", fontsize='small',fontname='Arial')
	    cbar2.ax.xaxis.set_label_position('top')
	#######################################
	resiTMP = medianRAresidual[cri,mjdIdx]
        ax = axs.ravel()[mjdIdx*3]
	img = ax.scatter(medianRADEC[cri,0], medianRADEC[cri,1], c=resiTMP, cmap=cmap0, norm=norm0, marker='o', s=15,lw=0, rasterized=True)
        #ax.set_xlim([(xlim[0]-xscale/20), (xlim[1]+xscale/20)])
        #ax.set_ylim([(ylim[0]-yscale/10), 25])#(ylim[1]+yscale/10)])
        ax.set_xlim([xmin, xmax])#([248.4, 252.5])
        ax.set_ylim([ymin, ymax])#([34.4, 38.5])
	ax.text(55, 1, r"\textbf{"+Years[mjdIdx]+"}", fontsize=16)
        ax.set_ylabel(r'\textbf{ $\delta$ (deg)}', fontsize=18)
        #ax.set_yticks([int(ylim[0]+yscale/10), int(ylim[0]+yscale/2), 22])#int(ylim[1]-yscale/10)])
	ax.set_yticks([0, 10 ,20])#([34, 36, 38]) #
	print np.median(medianMJD[cri,mjdIdx])
        print np.max(resiTMP), np.min(resiTMP)
        if((mjdIdx*3)>(mjdN*3-4)): 
	    ax.set_xlabel(r'\textbf{$\alpha$ (deg)}', fontsize=18)
	    ax.set_xticks([int(xlim[0]+xscale/10), int(xlim[0]+xscale/2), int(xlim[1]-xscale/10)])
	    #ax.set_xticks([248, 250, 252])
	################################
        ax = axs.ravel()[mjdIdx*3+1]
        stdTMP = stdRAresidual[cri,mjdIdx]
        img1 = ax.scatter(medianRADEC[cri,0], medianRADEC[cri,1], c=stdTMP, cmap=cmap1, norm=norm1, marker='o',s=15, lw=0, rasterized=True)
        ax.set_ylabel('')
        ax.set_yticks([])
        #ax.set_xlim([(xlim[0]-xscale/20), (xlim[1]+xscale/20)])
        #ax.set_ylim([(ylim[0]-yscale/10), 25])#(ylim[1]+yscale/10)])
        ax.set_xlim([xmin, xmax])#([248.4, 252.5])
        ax.set_ylim([ymin, ymax])#([34.4, 38.5])
        if((mjdIdx*3 + 1)>mjdN*3-4): 
	    ax.set_xlabel(r'\textbf{$\alpha$ (deg)}', fontsize=18)
	    ax.set_xticks([int(xlim[0]+xscale/10), int(xlim[0]+xscale/2), int(xlim[1]-xscale/10)])
	    #ax.set_xticks([248, 250, 252])
       #################################
	ax = axs.ravel()[mjdIdx*3+2]
        noTMP = noGalaxyMJD[cri,mjdIdx]
        img1 = ax.scatter(medianRADEC[cri,0], medianRADEC[cri,1], c=noTMP, cmap=cmap2, norm=norm2, marker='o',s=15, lw=0, rasterized=True)
	ax.plot(radec0.ra.value, radec0.dec.value, color="black", linewidth=2.5, linestyle="--")
	ax.plot(radec1.ra.value, radec1.dec.value, color="black", linewidth=2.5, linestyle="-")
	ax.plot(radec2.ra.value, radec2.dec.value, color="black", linewidth=2.5, linestyle="--")
        ax.set_ylabel('')
        ax.set_yticks([])
        #ax.set_xlim([(xlim[0]-xscale/20), (xlim[1]+xscale/20)])
        #ax.set_ylim([(ylim[0]-yscale/10), 25])#(ylim[1]+yscale/10)])
        ax.set_xlim([xmin, xmax])#([248.4, 252.5])
        ax.set_ylim([ymin, ymax])#([34.4, 38.5])
        if((mjdIdx*3 + 2)>mjdN*3-4): 
	    ax.set_xlabel(r'\textbf{$\alpha$ (deg)}', fontsize=18)
	    ax.set_xticks([int(xlim[0]+xscale/10), int(xlim[0]+xscale/2), int(xlim[1]-xscale/10)])
	    #ax.set_xticks([248, 250, 252])

    plt.savefig(root+"distri_galaxies_dec_PS1.pdf")



def pixelQSOs(fname, mr0, mr1, txt, root):
    qsopix = np.load(fname)
    x = qsopix['ra']-180
    y = qsopix['dec']

    degtorad = np.pi/180.
    xrad = x * degtorad
    yrad = y * degtorad
    if(txt[0]=="ra"):
    	sigpix = qsopix['epmraInPixel']
    elif(txt[0]=="dec"):
	sigpix = qsopix['epmdecInPixel']
    else:
	print "error in the txt[0] value, which only is 'ra' or 'dec' "
	return
    msigpix = np.median(sigpix)
    ssigpix = (0.741*(np.percentile(sigpix, 75) - np.percentile(sigpix, 25)))

    #############Galactic Disk###############
    lDisk = np.arange(0, 360, 1.0)
    bDisk0 = np.zeros(len(lDisk)) -20
    bDisk1 = np.zeros(len(lDisk)) + 0
    bDisk2 = np.zeros(len(lDisk)) + 20
    radec0 = SkyCoord("galactic", l=lDisk, b=bDisk0, unit=(u.degree, u.degree)).fk5
    radec1 = SkyCoord("galactic", l=lDisk, b=bDisk1, unit=(u.degree, u.degree)).fk5
    radec2 = SkyCoord("galactic", l=lDisk, b=bDisk2, unit=(u.degree, u.degree)).fk5


    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.2, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.16, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([msigpix-3*ssigpix, 0]), vmax=msigpix+3*ssigpix)#(vmin=0, vmax=5.0)#
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label(txt[3], fontsize=18,fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=sigpix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    ax1.set_xticklabels(np.arange(30,331,30))
    if(txt[0]=="ra"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"pmraerr.pdf")
    elif(txt[0]=="dec"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"pmdecerr.pdf")


    if(txt[0]=="ra"):
    	mupix = qsopix['mrapm']
    elif(txt[0]=="dec"):
	mupix = qsopix['mdecpm']
    else:
	print "error in the txt[0] value, which only is 'ra' or 'dec' "
	return
    mmupix = np.median(mupix)
    smupix = (0.741*(np.percentile(mupix, 75) - np.percentile(mupix, 25)))
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.2, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.16, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=mmupix - 3*smupix, vmax=mmupix + 3*smupix)#(vmin=-3.0, vmax=3.0)#
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8]) 
    cbar0.set_label(txt[4], fontsize=18,fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=mupix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    ax1.set_xticklabels(np.arange(30,331,30))
    if(txt[0]=="ra"):
        pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"mpmra.pdf")
    elif(txt[0]=="dec"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"mpmdec.pdf")

    #####another method to plot the projection, following code have some problem
    '''fig = pl.figure(figsize=(6.0,3.5))
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)
    phi = xrad + np.pi
    theta = np.pi/2 - yrad
    hp.visufunc.graticule()
    hp.visufunc.projscatter(theta, phi, lonlat=False, s = 20, coord='C', cmap = cmap0, norm=norm0, rot = (0, np.pi/2), c=mupix, marker='o', lw=0)
    pl.savefig(root+"figs/pixel"+label+"mpmra2.png")'''


    numpix = qsopix['numInRadius']
    mnumpix = np.median(numpix)
    snumpix = (0.741*(np.percentile(numpix, 75) - np.percentile(numpix, 25)))
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.2, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.16, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([mnumpix - 3*snumpix, 0]), vmax=mnumpix + 3*snumpix)
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label(r'$Number\ of\ '+txt[2]+'\ in\ each\ pixel$', fontsize=18,fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=numpix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    ax1.set_xticklabels(np.arange(30,331,30))
    pl.savefig(root+"figs/pixel"+txt[2]+"num.pdf")


def pixelQSOslb(fname, mr0, mr1, txt, root):
    qsopix = np.load(fname)
    ra = qsopix['ra']#-180
    dec = qsopix['dec']

    lb = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    l =lb.galactic.l.value
    b =lb.galactic.b.value

    degtorad = np.pi/180.
    xrad = (l-180) * degtorad
    yrad = b * degtorad
    if(txt[0]=="ra"):
    	sigpix = qsopix['epmraInPixel']
    elif(txt[0]=="dec"):
	sigpix = qsopix['epmdecInPixel']
    else:
	print "error in the txt[0] value, which only is 'ra' or 'dec' "
	return
    msigpix = np.median(sigpix)
    ssigpix = (0.741*(np.percentile(sigpix, 75) - np.percentile(sigpix, 25)))

    #############Galactic Disk###############
    lDisk = np.arange(0, 360, 1.0)
    bDisk0 = np.zeros(len(lDisk)) -20
    bDisk1 = np.zeros(len(lDisk)) + 0
    bDisk2 = np.zeros(len(lDisk)) + 20
    radec0 = bDisk0#SkyCoord("galactic", l=lDisk, b=bDisk0, unit=(u.degree, u.degree)).fk5
    radec1 = bDisk1#SkyCoord("galactic", l=lDisk, b=bDisk1, unit=(u.degree, u.degree)).fk5
    radec2 = bDisk2#SkyCoord("galactic", l=lDisk, b=bDisk2, unit=(u.degree, u.degree)).fk5


    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([msigpix-3*ssigpix, 0]), vmax=msigpix+3*ssigpix)#(vmin=0, vmax=5.0)#
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label(txt[3], fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=sigpix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    '''ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)'''
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    if(txt[0]=="ra"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"pmraerr_lb.pdf")
    elif(txt[0]=="dec"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"pmdecerr_lb.pdf")


    if(txt[0]=="ra"):
    	mupix = qsopix['mrapm']
    elif(txt[0]=="dec"):
	mupix = qsopix['mdecpm']
    else:
	print "error in the txt[0] value, which only is 'ra' or 'dec' "
	return
    mmupix = np.median(mupix)
    smupix = (0.741*(np.percentile(mupix, 75) - np.percentile(mupix, 25)))
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=mmupix - 3*smupix, vmax=mmupix + 3*smupix)#(vmin=-3.0, vmax=3.0)#
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8]) 
    cbar0.set_label(txt[4], fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=mupix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    '''ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)'''
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    if(txt[0]=="ra"):
        pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"mpmra_lb.pdf")
    elif(txt[0]=="dec"):
    	pl.savefig(root+"figs/pixel_"+txt[1]+"_"+txt[2]+"mpmdec_lb.pdf")

    #####another method to plot the projection, following code have some problem
    '''fig = pl.figure(figsize=(6.0,3.5))
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)
    phi = xrad + np.pi
    theta = np.pi/2 - yrad
    hp.visufunc.graticule()
    hp.visufunc.projscatter(theta, phi, lonlat=False, s = 20, coord='C', cmap = cmap0, norm=norm0, rot = (0, np.pi/2), c=mupix, marker='o', lw=0)
    pl.savefig(root+"figs/pixel"+label+"mpmra2.png")'''


    numpix = qsopix['numInRadius']
    mnumpix = np.median(numpix)
    snumpix = (0.741*(np.percentile(numpix, 75) - np.percentile(numpix, 25)))
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([mnumpix - 3*snumpix, 0]), vmax=mnumpix + 3*snumpix)
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label(r'$Number\ of\ '+txt[2]+'\ in\ each\ pixel$', fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=numpix, cmap=cmap0, norm=norm0, marker=u'o', s=20, lw=0)
    '''ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)'''
    ax1.grid(True)
    ax1.text(-60*degtorad, -55*degtorad, "$"+txt[1]+"\ "+txt[2]+"$", fontsize=10)
    ax1.text(30*degtorad, -55*degtorad, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=10)
    pl.savefig(root+"figs/pixel"+txt[2]+"num_lb.pdf")




def pixelOffSets(fname, ts, te, survey, root):
    qsopix = np.load(fname)
    x = qsopix['ra']-180
    y = qsopix['dec']

    degtorad = np.pi/180.
    xrad = x * degtorad
    yrad = y * degtorad
    sigpix = qsopix['eoffsetDecInPixel']
    #print sigpix

    #############Galactic Disk###############
    lDisk = np.arange(0, 360, 1.0)
    bDisk0 = np.zeros(len(lDisk)) -20
    bDisk1 = np.zeros(len(lDisk)) + 0
    bDisk2 = np.zeros(len(lDisk)) + 20
    radec0 = SkyCoord("galactic", l=lDisk, b=bDisk0, unit=(u.degree, u.degree)).fk5
    radec1 = SkyCoord("galactic", l=lDisk, b=bDisk1, unit=(u.degree, u.degree)).fk5
    radec2 = SkyCoord("galactic", l=lDisk, b=bDisk2, unit=(u.degree, u.degree)).fk5


    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    msigpix = np.median(sigpix)
    ssigpix = (0.741*(np.percentile(sigpix, 75) - np.percentile(sigpix, 25)))
    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([msigpix-3*ssigpix, 0]), vmax=msigpix+3*ssigpix)
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label('RMS of Offsets (DEC/mas)', fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=sigpix, cmap=cmap0, norm=norm0, marker='.', lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-100*degtorad, -55*degtorad, survey + " Galaxies", fontsize=10)
    ax1.text(10*degtorad, -55*degtorad, str(ts)[0:10]+" to "+str(te)[0:10], fontsize=10)
    pl.savefig(root+"figs/pixelOffSeterr" + survey +str(ts)[0:10]+str(te)[0:10]+".png")

    offsetpix = qsopix['moffsetDecInPixel']
    moffsetpix = np.median(offsetpix)
    soffsetpix = (0.741*(np.percentile(offsetpix-moffsetpix, 75) - np.percentile(offsetpix-moffsetpix, 25)))
    offsetpix = offsetpix-moffsetpix
    moffsetpix = np.median(offsetpix)
    #print mupix
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=moffsetpix-3*soffsetpix, vmax=moffsetpix+3*soffsetpix)
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label('Median Offsets (DEC/mas)', fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=offsetpix, cmap=cmap0, norm=norm0, marker='.', lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-100*degtorad, -55*degtorad, survey + " Galaxies", fontsize=10)
    ax1.text(10*degtorad, -55*degtorad, str(ts)[0:10]+" to "+str(te)[0:10], fontsize=10)
    pl.savefig(root+"figs/pixelmOffSet" + survey+str(ts)[0:10]+str(te)[0:10]+ ".png")

    numpix = qsopix['numInRadius']
    mnumpix = np.median(numpix)
    snumpix = (0.741*(np.percentile(numpix, 75) - np.percentile(numpix, 25)))
    #print numpix
    fig = pl.figure(figsize=(6.0,3.5))
    ax1 = fig.add_subplot(111, projection="mollweide")
    pl.subplots_adjust(left=0.08, bottom=0.18, right=0.96, top=0.99, wspace=0., hspace=0.)

    cax0 = fig.add_axes([0.05, 0.12, 0.9, 0.02])
    cmap0 = mpl.cm.jet  #cool
    norm0 = mpl.colors.Normalize(vmin=max([mnumpix-3*snumpix, 0]), vmax=mnumpix+3*snumpix)
    cbar0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap0,norm=norm0,orientation='horizontal')
    cax0.xaxis.set_ticks_position("bottom")
    #cbar0.set_ticks([-8,-4, 0,4, 8])
    cbar0.set_label('Galaxy Number in Each Pixel', fontsize='small',fontname='Arial')
    cbar0.ax.xaxis.set_label_position('bottom')
    ax1.scatter(xrad, yrad, c=numpix, cmap=cmap0, norm=norm0, marker='.', lw=0)
    ax1.plot((radec0.ra.value-180)*degtorad, radec0.dec.value*degtorad, "m,", ms=1.5)
    ax1.plot((radec1.ra.value-180)*degtorad, radec1.dec.value*degtorad, "m.", ms=1.5)
    ax1.plot((radec2.ra.value-180)*degtorad, radec2.dec.value*degtorad, "m,", ms=1.5)
    ax1.grid(True)
    ax1.text(-100*degtorad, -55*degtorad, survey + " Galaxies", fontsize=10)
    ax1.text(10*degtorad, -55*degtorad, str(ts)[0:10]+" to "+str(te)[0:10], fontsize=10)
    pl.savefig(root+"figs/pixelOffSetnum" + survey +str(ts)[0:10]+str(te)[0:10]+".png")
##############fitting PM in different objID################
def fitPMobj(obj_id, ra, raErr, mjd, mr, root):
    uniqueStar = np.unique(obj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    idx=0
    for idx in range(0, 1):#len(uniqueStar)):
	ind = obj_id == -4998010423962147576#uniqueStar[idx]
	tmp = ra[ind]
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = raErr[ind]*1000*3600
	xtmp = mjd[ind]
	print xtmp, ytmp, yErrtmp
	#xtmp = xtmp - xtmp[tmp == np.min(tmp)]	
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	#xtmp = xtmp[xtmp<56100]
	if((len(xtmp)>2)):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		#p,residuals, rank, singular_values, rcond = np.polyfit(xtmp, ytmp, 1, full=True)
		#x_extra = np.concatenate((xtmp,xtmp[-1:]))
		#y_extra = np.concatenate((ytmp, ytmp[-1:]))
		#weights = [1 for i in range(xtmp.size)]
		#weights.append(sys.float_info.epsilon)
		#p,cov = np.polyfit(xtmp, ytmp, 1, w=yErrtmp, cov=True)
		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		mu = p[0]*360
		muErr = np.sqrt(np.diag(cov))[0]*360
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		R2.append([objIDtmp, rMagtmp, Rsqr, mu, muErr])
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		ax.text(np.min(xtmp)+300,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
		ax.text(np.min(xtmp)+300,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)

		ax.plot(xtmp, xtmp*p[0]+p[1], '-')
		ax.set_ylim([int(-1*np.max(ytmp)/3), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_xlabel('MJD')
		ax.set_xlim(np.min(xtmp)-200,np.max(xtmp)+200) 
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.25, "objID = %d" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.15, "rmag = %0.2f" %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.05, "$\mu_{\eta}$ = %0.2f" %(mu)+"$\pm$%0.2f" \
		    %muErr, fontsize=10)
		ax.set_title('Postional Variation of Star in Different Epochs', fontsize=12)
		figname = "fig_PS1_Star" + "%d" %idx
		plt.savefig(root+"figs/%s_tst.png" %figname)
		#plt.show()

    R2 = np.array(R2,dtype=object)
    #np.save(root+"figs/R2_PS1_Gal_new2.npy", R2)


##############fitting PM in different objID################
def fitPMobjSDSS(obj_id, obj_idSDSS, ra, raErr, raSDSS, raErrSDSS, mjd, mjdSDSS, mr, root):
    print "len(obj_idSDSS):", len(obj_idSDSS)
    uniqueStar = np.unique(obj_idSDSS)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    idx=0
    for idx in range(0, len(uniqueStar)):
	ind = obj_id == uniqueStar[idx] #-6777811886076048565
	if (np.sum(ind)<2): continue
	indSDSS = obj_idSDSS == uniqueStar[idx]
	tmp = ra[ind]
	tmpSDSS = raSDSS[indSDSS]
	#print "tmpSDSS", tmpSDSS, tmp
	tmp = np.concatenate((tmp, tmpSDSS))
	ytmp = (tmp-np.min(tmp))*1000*3600
	#ytmp = (tmp)*1000*3600
	yErrtmp = np.concatenate((raErr[ind], raErrSDSS[indSDSS]))*1000*3600
	xtmp = mjd[ind]
	xtmp = np.concatenate((xtmp, mjdSDSS[indSDSS]))
	#print xtmp
	#xtmp = xtmp - xtmp[tmp == np.min(tmp)]
	
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	#xtmp = xtmp[xtmp<56100]
	if((len(xtmp)>2)):
		'''fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')'''
		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		mu = p[0]*360
		muErr = np.sqrt(np.diag(cov))[0]*360
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		#R2.append([rMagtmp, Rsqr, mu, muErr])
		#print "objIDtmp:", objIDtmp
		R2.append([objIDtmp, rMagtmp, Rsqr, mu, muErr])
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		'''ax.text(np.min(xtmp)+300,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
		ax.text(np.min(xtmp)+300,np.max(ytmp)*0.84,'$R^2 = %0.2f$'% Rsqr, fontsize=10)
		ax.plot(xtmp, xtmp*p[0]+p[1], '-')
		ax.set_ylim([int(-1*np.max(ytmp)/3), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_xlabel('$MJD')
		ax.set_xlim(np.min(xtmp)-200,np.max(xtmp)+200) 
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.15, "$rmag = %0.2f$" %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+300, np.max(ytmp)*1.05, "$\mu_{\eta} = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    %(np.sqrt(np.diag(cov))[0]*360), fontsize=10)
		ax.set_title('Postional Variation of Star in Different Epochs', fontsize=12)
		figname = "figSDSSStar_" + "%d" %idx
		plt.savefig(root+"figs/%s_25masErr.png" %figname)'''
		#plt.show()

    R2 = np.array(R2,dtype=object)# dtype=[('obj_id', 'i8'), ('rmag', 'f8'), ('R2', 'f8'), ('mu', 'f8'), ('muErr', 'f8')]
    np.save(root+"figs/R2_PS1_SDSS_GalETA.npy", R2)

##############fitting PM in different objID################
def fitPMobjPS1(obj_id, ra, raErr, mjd, mr, muPS1, muErrPS1, cobj_id, cra, craErr, cmjd, pflag, root):
    uniqueStar = np.unique(obj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    idx=0
    if(pflag):
	outliers = np.array([-4799869632544364352, \
	    -4763770466781189284, -4997992831776515880])# -4799869632544364323,-4763770466781242308, 
	objID_curr = outliers
	N = len(outliers)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	ind = obj_id == objID_curr[idx] #-6777811886076048565 #
	xtmp, uniqueIdx = np.unique(mjd[ind], return_index=True)
	print xtmp
	#print "obj_id[ind]:", obj_id[ind]
	#xtmp = (mjd[ind])
	tmp = ra[ind]
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = raErr[ind]*1000*3600
	ytmp = ytmp[uniqueIdx]
	yErrtmp = yErrtmp[uniqueIdx]
	#print len(xtmp), len(ytmp), len(yErrtmp)
	#xtmp = xtmp - xtmp[tmp == np.min(tmp)]	
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	muPS1tmp = muPS1[ind][0]
	muErrPS1tmp = muErrPS1[ind][0]

	cind = cobj_id == objID_curr[idx]#-4997992831776515880#
	cxtmp = (cmjd[cind])
	ctmp = cra[cind]
	#print "min(tmp)-(ctmp)=", (np.min(tmp)-np.min(ctmp))*1000*3600
	cytmp = (ctmp-np.min(tmp))*1000*3600
	cyErrtmp = craErr[cind]*1000*3600

	if((len(xtmp)>2) and (len(cxtmp)>2)):		
	        #print cxtmp, cytmp, cyErrtmp
	        cp, ccov = curve_fit(func, cxtmp, cytmp, sigma=cyErrtmp, absolute_sigma=True)
	        cmu = cp[0]*360
	        cmuErr = np.sqrt(np.diag(ccov))[0]*360

		if(pflag):
			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.errorbar(cxtmp, cytmp, cyErrtmp, fmt='.', color='g')
			ax.plot(xtmp, xtmp*cp[0]+cp[1], 'g-')
			ax.plot(xtmp, xtmp*(muPS1tmp/360.0)+(np.median(ytmp)-muPS1tmp/360.0*np.median(xtmp)), 'k--')

		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		mu = p[0]*360
		muErr = np.sqrt(np.diag(cov))[0]*360
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		R2.append([objIDtmp, rMagtmp, Rsqr, mu, muErr, muPS1tmp, muErrPS1tmp, cmu, cmuErr])
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		if(pflag):
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)

			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel('$\eta-min(\eta)/mas$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ", $\mu_{\delta}(PS1) = %0.2f$" %(muPS1tmp)+"$\pm%0.2f$" %muErrPS1tmp \
		    	    + ", $\mu_{\delta}(AVG) = %0.2f$" %(cmu)+"$\pm%0.2f$" %cmuErr, fontsize=10)
			ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
			figname = "fig_PS1_original_combined" + "%d" %idx
			plt.savefig(root+"figs/%s_outlier_cpm_h.png" %figname)
			#plt.show()
    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/R2_PS1_original_new.npy", R2)

##############fitting PM in different objID################
def fitPMobjPS12(obj_id, ra, raErr, mjd, mr, cpm, cobj_id, cra, craErr, cmjd, pflag, root):
    uniqueStar = np.unique(obj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    idx=0
    if(pflag):
	outliers = np.array([-4799869632544364352, -4997992831776515880, -4763770466781189284, -4763770466781242308, \
	    -4799869632544364323])#
	objID_curr = outliers
	N = len(outliers)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	ind = obj_id == objID_curr[idx] #-6777811886076048565 #
	xtmp, uniqueIdx = np.unique(mjd[ind], return_index=True)
	print xtmp
	#print "obj_id[ind]:", obj_id[ind]
	#xtmp = (mjd[ind])
	tmp = ra[ind]
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = raErr[ind]*1000*3600
	ytmp = ytmp[uniqueIdx]
	yErrtmp = yErrtmp[uniqueIdx]
	#print len(xtmp), len(ytmp), len(yErrtmp)
	#xtmp = xtmp - xtmp[tmp == np.min(tmp)]	
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	ind4PS1 = cpm['obj_id']==objID_curr[idx]
	muPS1tmp = cpm['U_DEC'][ind4PS1][0]*1000
	muErrPS1tmp = cpm['V_DEC_ERR'][ind4PS1][0]*1000

	cind = cobj_id == objID_curr[idx]#-4997992831776515880#
	cxtmp = (cmjd[cind])
	ctmp = cra[cind]
	#print "min(tmp)-(ctmp)=", (np.min(tmp)-np.min(ctmp))*1000*3600
	cytmp = (ctmp-np.min(tmp))*1000*3600
	cyErrtmp = craErr[cind]*1000*3600

	if((len(xtmp)>2) and (len(cxtmp)>2)):		
	        #print cxtmp, cytmp, cyErrtmp
	        cp, ccov = curve_fit(func, cxtmp, cytmp, sigma=cyErrtmp, absolute_sigma=True)
	        cmu = cp[0]*360
	        cmuErr = np.sqrt(np.diag(ccov))[0]*360

		if(pflag):
			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.errorbar(cxtmp, cytmp, cyErrtmp, fmt='.', color='g')
			ax.plot(xtmp, xtmp*cp[0]+cp[1], 'g-')
			ax.plot(xtmp, xtmp*(muPS1tmp/360.0)+(np.median(ytmp)-muPS1tmp/360.0*np.median(xtmp)), 'k--')

		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		mu = p[0]*360
		muErr = np.sqrt(np.diag(cov))[0]*360
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		R2.append([objIDtmp, rMagtmp, Rsqr, mu, muErr, muPS1tmp, muErrPS1tmp, cmu, cmuErr])
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		if(pflag):
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)

			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel('$\eta-min(\eta)/mas$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ", $\mu_{\delta}(PS1) = %0.2f$" %(muPS1tmp)+"$\pm%0.2f$" %muErrPS1tmp \
		    	    + ", $\mu_{\delta}(AVG) = %0.2f$" %(cmu)+"$\pm%0.2f$" %cmuErr, fontsize=10)
			ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
			figname = "fig_PS1_original_combined" + "%d" %idx
			plt.savefig(root+"figs/%s_outlier_cpm2_no_offsetErr.png" %figname)
			#plt.show()
    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/R2_PS1_original_new.npy", R2)

##############fitting PM in different objID################
def fitPMobj2MASS(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, raErr2, mjd2, obj_id3, ra3, raErr3, mjd3, cobj_id, cmuDEC, cmuDECerr, pflag, root):
    uniqueStar = np.unique(obj_id0)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    if(pflag):
	outliers = np.array([-4997992831776515880, -4799869632544364352, -4763770466781189284, -4763770466781242308, \
	    -4799869632544364323])#
	objID_curr = uniqueStar#outliers
	N = len(uniqueStar)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	 
	ind1 = obj_id1 == objID_curr[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	

	ind0 = obj_id0 == objID_curr[idx]
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]

	ind3 = obj_id3 == objID_curr[idx]
	xtmp3, uniqueIdx3 = np.unique(mjd3[ind3], return_index=True)
	tmp3 = ra3[ind3][uniqueIdx3]

	tmp = np.concatenate((tmp0, tmp1, tmp3))
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], raErr1[ind1][uniqueIdx1], raErr3[ind3][uniqueIdx3]))*1000*3600
	xtmp = np.concatenate((xtmp0, xtmp1, xtmp3))
	#print 'xtmp, ytmp,yErrtmp', xtmp,ytmp,yErrtmp

	objIDtmp = obj_id0[ind0][0]
	rMagtmp = mr0[ind0][0]
	ind4PS1 = cobj_id==objID_curr[idx]
	if np.sum(ind4PS1)>0: 
	    muPS1tmp = cmuDEC[ind4PS1][0]*1000
	    muErrPS1tmp = cmuDECerr[ind4PS1][0]*1000
	else:
	    muPS1tmp = 999
	    muErrPS1tmp = 999

	ind2 = obj_id2 == objID_curr[idx]
	xtmp2 = (mjd2[ind2])
	xtmp2 = np.concatenate((xtmp0, xtmp2, xtmp3))

	tmp2 = ra2[ind2]
	tmp2 = np.concatenate((tmp0, tmp2, tmp3))

	ytmp2 = (tmp2-np.min(tmp))*1000*3600
	yErrtmp2 = np.concatenate((raErr0[ind0][uniqueIdx0], raErr2[ind2], raErr3[ind3][uniqueIdx3]))*1000*3600
	
	# count observational times in each season for one star
	nObsTmp = np.zeros(len(mjd2[ind2])) 
	if np.sum(ind1)>0: 
	    for mjd1Idx in range(0, len(mjd2[ind2])):
	        mjdIndex = (xtmp1>mjd2[ind2][mjd1Idx]-70) & (xtmp1<mjd2[ind2][mjd1Idx]+70)
	        nObsTmp[mjd1Idx] = np.sum(mjdIndex)
	print "nObsTmp:", nObsTmp
	if((len(xtmp)>2) and (len(xtmp2)>2)):		
	        #print cxtmp, cytmp, cyErrtmp
	        p2, cov2 = curve_fit(func, xtmp2, ytmp2, sigma=yErrtmp2, absolute_sigma=True)
	        mu2 = p2[0]*365
	        muErr2 = np.sqrt(np.diag(cov2))[0]*365

		num2 = len(raErr2[ind2]) # Here num2 means the num of PS1
		num0 = len(raErr0[ind0][uniqueIdx0]) # Here num0 means the num of 2MASS	
		print "num0:", num0	
		dyX = np.zeros(num2)
		for idxX in range(0, num2): # X means cross-validation
			yErrtmp2[num0+idxX] = yErrtmp2[num0+idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
			p2X, cov2X = curve_fit(func, xtmp2, ytmp2, sigma=yErrtmp2, absolute_sigma=True)
			dyX[idxX] = (ytmp2[num0+idxX] - (p2X[0]*xtmp2[num0+idxX] + p2X[1]))
			yErrtmp2[num0+idxX] = yErrtmp2[num0+idxX] - 9999.0 #Recover the current uncertainty
		#print "dyX:", dyX	
		if(pflag):
			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.errorbar(xtmp2, ytmp2, yErrtmp2, fmt='.', color='g')
			ax.plot(xtmp, xtmp*p2[0]+p2[1], 'g-')
			ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')

		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		mu = p[0]*365
		muErr = np.sqrt(np.diag(cov))[0]*365
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		R2.append([objIDtmp, rMagtmp, Rsqr, mu, muErr, muPS1tmp, muErrPS1tmp, mu2, muErr2, dyX, yErrtmp2[num0:(num0+num2)], nObsTmp])
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		if(pflag):
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)

			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel('$\eta-min(\eta)/mas$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ", $\mu_{\delta}(PS1) = %0.2f$" %(muPS1tmp)+"$\pm%0.2f$" %muErrPS1tmp \
		    	    + ", $\mu_{\delta}(AVG) = %0.2f$" %(mu2)+"$\pm%0.2f$" %muErr2, fontsize=10)
			ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
			figname = "fig_PS1_2MASS_SDSS" + "%d" %idx
			plt.savefig(root+"figs/ra13p613p8_dec34p6534p85/%s.png" %figname)
			#plt.show()
    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/R2_PS1_2MASS_PS1_X.npy", R2)

##############fitting PM in different objID, , This is for GPS1 case################
def fitPMobjPal5(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, \
	raErr2, mjd2, pobj_id, pmuDEC, pmuDECerr, tobj_id, tmuDEC, tmuDECerr,  \
	gobj_id, gra, graErr, gmjd, cpm, pflag, root):
    uniqueStar = np.unique(tobj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    if(pflag):
	outliers = np.array([4796122496917232244, 4796157681289161166, 4796104904731190081, 4796122496917242833, 4796157681289159480, 4796157681289159647, 4796157681289376695, 4814154487612597084]) #
	objID_curr = outliers#
	N = 1#len(outliers)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	 
	ind1 = obj_id1 == objID_curr[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	#print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == objID_curr[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]
	print "obj_id of PS1:", obj_id0[ind0][0], tmp0

	ind2 = obj_id2 == objID_curr[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]
	print xtmp2, xtmp1, tmp2, tmp1

	ind3 = cpm['obj_id'] == objID_curr[idx]
	xtmp3, uniqueIdx3 = np.unique(cpm['mjd'][ind3], return_index=True)
	tmp3 = cpm['dec'][ind3][uniqueIdx3]
	yErrtmp3 = cpm['decErr'][ind3][uniqueIdx3]*1000*3600
	band = cpm['band'][ind3][uniqueIdx3]
	x_fix =  cpm['x_fix'][ind3][uniqueIdx3]
	y_fix =  cpm['y_fix'][ind3][uniqueIdx3]
	#print xtmp3, tmp3

	gind = gobj_id == objID_curr[idx]
	if (np.sum(gind)<1): continue
	gxtmp, guniqueIdx = np.unique(gmjd[gind], return_index=True)
	gtmp = gra[gind][guniqueIdx]
	print "obj_id of GAIA:", gobj_id[gind]

	tmp = np.concatenate((tmp1, tmp2, tmp0, gtmp))#tmp1, 
	#tmp = tmp0
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = np.concatenate(( raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2], raErr0[ind0][uniqueIdx0], graErr[gind][guniqueIdx]))*1000*3600
	#yErrtmp = raErr0[ind0][uniqueIdx0]*1000*3600
	xtmp = np.concatenate((xtmp1, xtmp2, xtmp0, gxtmp))# 
	#xtmp = xtmp0
	#print 'xtmp, ytmp,yErrtmp', xtmp,ytmp,yErrtmp

	tmpCMP = np.concatenate((tmp1, tmp2, tmp3, gtmp))#tmp1, 
	ytmpCMP = (tmpCMP-np.min(tmp))*1000*3600
	yErrtmpCMP = np.concatenate((raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2], cpm['decErr'][ind3][uniqueIdx3], graErr[gind][guniqueIdx]))*1000*3600 
	#yErrtmp = raErr0[ind0][uniqueIdx0]*1000*3600
	xtmpCMP = np.concatenate((xtmp1, xtmp2, xtmp3, gxtmp))

	tmp3 = (tmp3 - np.min(tmp))*1000*3600


	ind4PS1 = pobj_id==objID_curr[idx]
	if np.sum(ind4PS1)>0: 
	    muPS1tmp = pmuDEC[ind4PS1][0]*1000
	    muErrPS1tmp = pmuDECerr[ind4PS1][0]*1000
	else:
	    muPS1tmp = 999
	    muErrPS1tmp = 999

	ind4Pal5 = tobj_id==objID_curr[idx]
	if np.sum(ind4Pal5)>0: 
	    muPal5tmp = tmuDEC[ind4Pal5][0]
	    muErrPal5tmp = tmuDECerr[ind4Pal5][0]
	else:
	    muPal5tmp = 999
	    muErrPal5tmp = 999
	
	if((len(xtmp)>2)):
	    objIDtmp = objID_curr[idx]
	    rMagtmp = mr0[ind0][0]
	    s0 = len(xtmp1) + len(xtmp2)
	    num0 = len(raErr0[ind0][uniqueIdx0]) + 1
	    xpm = np.zeros([5, num0])
	    for idxX in range(0, num0): # X means cross-validation
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] + 99999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = abs(ytmp[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - np.delete(ytmp, s0+idxX))**2/np.delete(yErrtmp, s0+idxX)**2)/(xtmp.size-2)
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] - 99999.0


	    gp, gcov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    gmu = gp[0]*365 
	    gmuErr = np.sqrt(np.diag(gcov))[0]*365
	    gchi2 = np.sum((gp[0]*xtmp + gp[1] - ytmp)**2/yErrtmp**2)/(xtmp.size-1)
	    chi2Tmp = np.concatenate((xpm[4], np.array([gchi2])))
	    gmuTmp = np.concatenate((xpm[0], np.array([gmu])))
	    gslpTmp = np.concatenate((xpm[3], np.array([gp[1]])))
 	    idxO = np.argmin(chi2Tmp) 
	    gmuX = gmuTmp[idxO]
	    gslp = gslpTmp[idxO]
	    gchi2 = np.min(chi2Tmp)


	    mu = xpm[0][num0-1]  # No GAIA point case, which need GAIA point in the last
	    muErr = xpm[1][num0-1]
	    chi2 = xpm[4][num0-1]
			
	    ddy = (tmp3 - (gp[0]*xtmp3 + gp[1]))
	    offset = [ddy, band, x_fix, y_fix]

	    pCMP, covCMP = curve_fit(func, xtmpCMP, ytmpCMP, sigma=yErrtmpCMP, absolute_sigma=True)
	    muCMP = pCMP[0]*365
	    muErrCMP = np.sqrt(np.diag(covCMP))[0]*365
	    chi2CMP = np.sum((pCMP[0]*xtmpCMP + pCMP[1] - ytmpCMP)**2/yErrtmpCMP**2)/(xtmpCMP.size-1)

	    #R2.append([objIDtmp, rMagtmp, gchi2, gmu, gmuErr, muPS1tmp, muErrPS1tmp, muPal5tmp, muErrPal5tmp, xpm, offset])


	    if(pflag):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp3, tmp3, yErrtmp3, fmt='.', color='blue')
		ax.errorbar(xtmp[s0:s0+num0], ytmp[s0:s0+num0], yErrtmp[s0:s0+num0], fmt='o', color='red', elinewidth=2)
		ax.errorbar(xtmp2, (tmp2-min(tmp))*3600000.0, raErr2[ind2][uniqueIdx2]*3600000, fmt='o', color='m', elinewidth=2)
		ax.errorbar(gxtmp, (gtmp-min(tmp))*3600000.0, graErr[gind][guniqueIdx]*3600000, fmt='o', color='y', elinewidth=2)
		ax.errorbar(xtmp1, (tmp1-min(tmp))*3600000.0, raErr1[ind1][uniqueIdx1]*3600000, fmt='o', color='k', elinewidth=2)
		ax.plot(xtmp, xtmp*(muPal5tmp/365.0)+(np.median(ytmp)-muPal5tmp/365.0*np.median(xtmp)), 'g--')
		ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')

		ax.text(np.min(xtmp)+550,np.max(ytmp)*1.0,'$\chi^2 = %0.2f$'% gchi2 + '$(%0.2f$)'% chi2CMP, fontsize=10)
		#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)
		ax.plot(xtmpCMP, xtmpCMP*pCMP[0]+pCMP[1], 'b-')
		ax.plot(xtmp, xtmp*gmuX/365.0+gslp, 'r-')
		ax.set_ylim([int(-1*np.max(ytmp)/4)-10, np.max(ytmp)*1.4+5])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		#ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_ylabel('$\delta-min(\delta)/mas$')
		ax.set_xlabel('$MJD$')
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

		ax.text(np.min(xtmp)+550, np.max(ytmp)*1.30, "$objID = %d$" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+550, np.max(ytmp)*1.20, "$rmag = %0.2f$" \
		    %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+550, np.max(ytmp)*1.1, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" \
		    %gmuErr + "$(%0.2f$" %(muCMP)+"$\pm%0.2f)$" %muErrCMP \
		    + ", $\mu_{\delta, o} = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp + \
		    "$(%0.2f$" %(muPS1tmp)+"$\pm%0.2f)$" %muErrPS1tmp, fontsize=10)
		ax.set_title('Comparision of Proper Motions With Different Fitting Methods', fontsize=12)
		ax.minorticks_on()
		figname = "fig_PS1_Pal5_CPM_" + "%d" %idx
		plt.savefig(root+"figs/%s_GPS1.eps" %figname)

    #if(not pflag):
    	#R2 = np.array(R2,dtype=object)
    	#np.save(root+"figs/R2_OnlyPS1_CPM_Pal5_ppmxl.npy", R2)


##############fitting PM in different objID, This is for GP1 case################
def fitPMobjPal5b(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, \
	raErr2, mjd2, pobj_id, pmuDEC, pmuDECerr, tobj_id, tmuDEC, tmuDECerr,  \
	gobj_id, gra, graErr, gmjd, cpm, GP, pflag, root):
    uniqueStar = np.unique(tobj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    if(pflag):
	outliers = np.array([4796122496917232282, 4796122496917232287, 4796122496917232261, \
	    4796104904731181545, 4796122496917232244, 4796157681289161166, 4814154487612597084]) #
	objID_curr = outliers#uniqueStar#
	N = 3#len(outliers)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	 
	ind1 = obj_id1 == objID_curr[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	#print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == objID_curr[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]
	print "obj_id of PS1:", obj_id0[ind0][0], tmp0

	ind2 = obj_id2 == objID_curr[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]
	print xtmp2, xtmp1, tmp2, tmp1

	ind3 = cpm['obj_id'] == objID_curr[idx]
	xtmp3, uniqueIdx3 = np.unique(cpm['mjd'][ind3], return_index=True)
	tmp3 = cpm['dec'][ind3][uniqueIdx3]
	yErrtmp3 = cpm['decErr'][ind3][uniqueIdx3]*1000*3600
	band = cpm['band'][ind3][uniqueIdx3]
	x_fix =  cpm['x_fix'][ind3][uniqueIdx3]
	y_fix =  cpm['y_fix'][ind3][uniqueIdx3]
	#print xtmp3, tmp3

	gind = gobj_id == objID_curr[idx]
	if (np.sum(gind)<1): continue
	gxtmp, guniqueIdx = np.unique(gmjd[gind], return_index=True)
	gtmp = gra[gind][guniqueIdx]
	print "obj_id of GAIA:", gobj_id[gind]

	if(GP):
	    tmp = np.concatenate((tmp0, gtmp))# 
	    yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], graErr[gind][guniqueIdx]))*1000*3600
	    xtmp = np.concatenate((xtmp0, gxtmp))
	    tmpCMP = np.concatenate((tmp3, gtmp))
	    yErrtmpCMP = np.concatenate((cpm['decErr'][ind3][uniqueIdx3], graErr[gind][guniqueIdx]))*1000*3600
	    xtmpCMP = np.concatenate((xtmp3, gxtmp))
	else:
	    tmp = np.concatenate((tmp2,tmp1,tmp0, gtmp))
	    yErrtmp = np.concatenate((raErr2[ind2][uniqueIdx2], raErr1[ind1][uniqueIdx1], \
		raErr0[ind0][uniqueIdx0], graErr[gind][guniqueIdx]))*1000*3600	
	    xtmp = np.concatenate((xtmp2,xtmp1, xtmp0, gxtmp))
	    tmpCMP = np.concatenate((tmp2,tmp1,tmp3, gtmp))
	    yErrtmpCMP = np.concatenate((raErr2[ind2][uniqueIdx2], raErr1[ind1][uniqueIdx1], \
		cpm['decErr'][ind3][uniqueIdx3], graErr[gind][guniqueIdx]))*1000*3600
	    xtmpCMP = np.concatenate((xtmp2, xtmp1, xtmp3, gxtmp))
	
	ytmp = (tmp-np.min(tmp))*1000*3600 
	ytmpCMP = (tmpCMP-np.min(tmp))*1000*3600
	#tmp3 = (tmp3 - np.min(tmp))*1000*3600
	#tmp0 = (tmp0 - np.min(tmp))*1000*3600



	ind4PS1 = pobj_id==objID_curr[idx]
	if np.sum(ind4PS1)>0: 
	    muPS1tmp = pmuDEC[ind4PS1][0]*1000
	    muErrPS1tmp = pmuDECerr[ind4PS1][0]*1000
	else:
	    muPS1tmp = 999
	    muErrPS1tmp = 999

	ind4Pal5 = tobj_id==objID_curr[idx]
	if np.sum(ind4Pal5)>0: 
	    muPal5tmp = tmuDEC[ind4Pal5][0]
	    muErrPal5tmp = tmuDECerr[ind4Pal5][0]
	else:
	    muPal5tmp = 999
	    muErrPal5tmp = 999
	
	if((len(xtmp)>2)):
	    objIDtmp = objID_curr[idx]
	    rMagtmp = mr0[ind0][0]
	    s0 = 0
	    if(GP): num0 = len(xtmp0)  + len(gxtmp) #X-validate PS1 + GAIA + SDSS + 2MASS points
	    else: num0 = len(xtmp0)  + len(gxtmp) + len(xtmp1) + len(xtmp2)
	    xpm = np.zeros([5, num0])
	    for idxX in range(0, num0): # X means cross-validation
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] + 99999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = abs(ytmp[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - np.delete(ytmp, s0+idxX))**2/np.delete(yErrtmp, s0+idxX)**2)/(xtmp.size-3)
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] - 99999.0


	    gp, gcov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    gmu = gp[0]*365 
	    gmuErr = np.sqrt(np.diag(gcov))[0]*365
	    gchi2 = np.sum((gp[0]*xtmp + gp[1] - ytmp)**2/yErrtmp**2)/(xtmp.size-2)
	    chi2Tmp = np.concatenate((xpm[4], np.array([gchi2])))
	    gmuTmp = np.concatenate((xpm[0], np.array([gmu])))
	    gslpTmp = np.concatenate((xpm[3], np.array([gp[1]])))
 	    idxO = np.argmin(chi2Tmp) 
	    gmuX = gmuTmp[idxO]
	    gslp = gslpTmp[idxO]
	    gchi2 = np.min(chi2Tmp)


	    mu = xpm[0][num0-1]  # No GAIA point case, which need GAIA point in the last
	    muErr = xpm[1][num0-1]
	    slp = xpm[3][num0-1]
	    chi2 = xpm[4][num0-1]
			

	    pCMP, covCMP = curve_fit(func, xtmpCMP, ytmpCMP, sigma=yErrtmpCMP, absolute_sigma=True)
	    muCMP = pCMP[0]*365
	    muErrCMP = np.sqrt(np.diag(covCMP))[0]*365
	    chi2CMP = np.sum((pCMP[0]*xtmpCMP + pCMP[1] - ytmpCMP)**2/yErrtmpCMP**2)/(xtmpCMP.size-2)


	    if(pflag):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.15, bottom=0.15, right=0.96, top=0.91, wspace=0., hspace=0.)
		if(not GP): ax.errorbar(xtmp3, (tmp3-min(tmp))*3600000.0, yErrtmp3, fmt='.', color='blue')
		ax.errorbar(xtmp0, (tmp0-min(tmp))*3600000.0, raErr0[ind0][uniqueIdx0]*3600000, fmt='o', color='red', elinewidth=2)
		if(not GP): ax.errorbar(xtmp2, (tmp2-min(tmp))*3600000.0, raErr2[ind2][uniqueIdx2]*3600000, fmt='o', color='m', elinewidth=2)
		ax.errorbar(gxtmp, (gtmp-min(tmp))*3600000.0, graErr[gind][guniqueIdx]*3600000, fmt='o', color='y', elinewidth=2)
		if(not GP): ax.errorbar(xtmp1, (tmp1-min(tmp))*3600000.0, raErr1[ind1][uniqueIdx1]*3600000, fmt='o', color='k', elinewidth=2)
		ax.plot(xtmp, xtmp*(muPal5tmp/365.0)+(np.median(ytmp)-muPal5tmp/365.0*np.median(xtmp)), 'g--')
		if(not GP): ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')

		
		#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% Rsqr, fontsize=10)
		#ax.plot(xtmp, xtmp*gp[0]+gp[1], 'r--') # with gaia, but no removing outlier
		if(GP): ax.plot(xtmp, xtmp*mu/365.0+slp, 'r--') # no gaia, no removing outlier
		ax.plot(xtmp, xtmp*gmuX/365.0+gslp, 'r-')
		if(not GP): ax.plot(xtmpCMP, xtmpCMP*pCMP[0]+pCMP[1], 'b-')
		ax.set_ylim([int(-1*np.max(ytmp)/3), np.max(ytmp)*1.5])	
		ax.set_yticks([0, int(np.max(ytmp)/2)-1, int(np.max(ytmp))-2])
		#ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_ylabel('$(\delta-min(\delta))/mas$', fontsize=18)
		ax.set_xlabel('$MJD$', fontsize=18)
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

		if(GP): 
		    ax.text(np.min(xtmp)+100, np.max(ytmp)*1.40, "$objID = %d$" %objIDtmp, fontsize=10)
		    ax.text(np.min(xtmp)+100, np.max(ytmp)*1.30, "$rmag = %0.2f$" \
		        %np.round(rMagtmp,decimals=2), fontsize=10)
		    ax.text(np.min(xtmp)+100,np.max(ytmp)*1.1,'$\chi_{\\nu}^2 = %0.2f$'% gchi2 + '$(%0.2f$)'% chi2, fontsize=10) #chi2Tmp[-1]
		    ax.text(np.min(xtmp)+100, np.max(ytmp)*1.2, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" \
		        %gmuErr + "$(%0.2f$" %(mu)+"$\pm%0.2f)$" %muErr \
		        + ", $\mu_{\delta, o} = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp \
		        , fontsize=10)
		else:
		    ax.text(np.min(xtmp)+550, np.max(ytmp)*1.40, "$objID = %d$" %objIDtmp, fontsize=10)
		    ax.text(np.min(xtmp)+550, np.max(ytmp)*1.30, "$rmag = %0.2f$" \
		        %np.round(rMagtmp,decimals=2), fontsize=10)
		    ax.text(np.min(xtmp)+550,np.max(ytmp)*1.1,'$\chi_{\\nu}^2 = %0.2f$'% gchi2 + '$(%0.2f$)'% chi2CMP, fontsize=10)
		    ax.text(np.min(xtmp)+550, np.max(ytmp)*1.2, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" \
		        %gmuErr + "$(%0.2f$" %(muCMP)+"$\pm%0.2f)$" %muErrCMP \
		        + ", $\mu_{\delta, o} = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp + \
		        "$(%0.2f$" %(muPS1tmp)+"$\pm%0.2f)$" %muErrPS1tmp, fontsize=10)
		ax.set_title('Comparision of Proper Motions With Different Fit Methods', fontsize=12)
		ax.minorticks_on()
		if(GP): figname = "fig_PS1_Pal5_CPM_" + "%d" %idx
		else: figname = "fig_GPS1_Pal5_CPM_" + "%d" %idx
		plt.savefig(root+"figs/%s_mrcal.eps" %figname)
##############fitting PM in different objID################
def fitPMobjGAIA(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, raErr2, mjd2, \
    tobj_id, tmuRA, tmuRAerr, tra, tdec, pflag, root):

    ################ create offset for each GAIA star #############
    nside = 2**6 # for GAIA offset, pixel area: 0.83929364521116678 deg^2 
    phiForObj   = (tra*np.pi)/180 
    thetaForObj = (90 - tdec)* (np.pi/180) 
    pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 	    
    pixelIndexArray = np.unique(pixelIndexForObj)
    Npixel = pixelIndexArray.size
    print "Going to process %d pixels." % Npixel
    offsetArray = np.zeros(len(tra))
    offsets = np.random.random(Npixel)*20 - 10 #[-10,10]
    for idx in range(0, Npixel):
        indexInPixel = pixelIndexForObj == pixelIndexArray[idx]
        offsetArray[indexInPixel] = offsets[idx] # it is coresponding to the tojb_id one by one


    gobj_id = np.zeros(len(tra), dtype='i8')
    gra = np.zeros(len(tra))
    #gdec = np.zeros(len(tra))
    gmjd = np.zeros(len(tra))
    ################ asign the position for each GAIA star #############
    N=len(tobj_id)
    #objID_curr = tobj_id
    for idx in range(0, N):
	 
	ind1 = obj_id1 == tobj_id[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	#print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == tobj_id[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]

	ind2 = obj_id2 == tobj_id[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]
	#print xtmp2, xtmp1, tmp2, tmp1

	tmp = np.concatenate((tmp0, tmp1, tmp2))#
	#tmp = tmp0
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2]))*1000*3600 #, 
	#yErrtmp = raErr0[ind0][uniqueIdx0]*1000*3600
	xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))#
	
	if((len(xtmp)>2)):
	    num0 = len(raErr0[ind0][uniqueIdx0]) # Here num0 means the num of PS1	
	    print "num0:", num0	
	    xpm = np.zeros([4, num0])
	    for idxX in range(0, num0): # X means cross-validation
		yErrtmp[idxX] = yErrtmp[idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = (ytmp[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		yErrtmp[idxX] = yErrtmp[idxX] - 9999.0

	    #print "muX", idx, ":", np.median(xpm[0])
	    #idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
	    #muX = np.median(np.delete(xpm[0], np.argmax(abs(xpm[2])))) # (Here should be wrong, since removing the most import value)
	    sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
	    idxO = abs(xpm[2] - np.median(xpm[2]))>1.5*sigTmp
	    if(sum(idxO)):
		idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
		gmuX = xpm[0][idxO]
	        gslp = xpm[3][idxO]
	    else:
		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		gmuX = p[0]*365.0
		gslp = p[1]
		#gmuX = np.median(xpm[0]) # This case does not work very well
		#gslp = (np.median(ytmp)-gmuX/365.0*np.median(xtmp))
	    #yErrtmp[idxO] = yErrtmp[idxO] + 9999.0
	    #p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    #yErrtmp[idxO] = yErrtmp[idxO] - 9999.0
	    gmjd[idx] = 57174  # 2015-06-01
	    gobj_id[idx] = tobj_id[idx].astype('i8')
	    #print gmjd[idx], muX, slp, offsetArray[idx]
	    gra[idx] = (gmjd[idx] *gmuX/365.0 + gslp + offsetArray[idx])/(1000*3600.0) + np.min(tmp) # predict GAIA's positon and adding an offset (deg)

    gra = gra + np.random.normal(0,3.0,N)/3600000.0 # adding 3mas observational error for each star



    ################# recovering GAIA's offset #############
    nside = 2**8 # pixel area: 0.0032784908016061202 deg^2
    phiForObj   = (tra*np.pi)/180 
    thetaForObj = (90 - tdec)* (np.pi/180) 
    pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 	    
    pixelIndexArray = np.unique(pixelIndexForObj)
    Npixel = pixelIndexArray.size
    print "Going to process %d pixels." % Npixel
    theta, phi = hp.pix2ang(nside, pixelIndexArray)
    pixelRa  = 180*phi/np.pi  # ra and dec of the centers of cells
    pixelDec = 90 - theta*180/np.pi
    cN = 20

    roffsetArray = np.zeros(len(tra))  # recovered offset
    for idx in range(0, Npixel):
        distTmp = sphdist(pixelRa[idx], pixelDec[idx], tra, tdec)
        angSepMask = np.argsort(distTmp)
        angSepMask = angSepMask[0:cN]
	objInRadius = tobj_id[angSepMask]
	indInPixel = pixelIndexForObj == pixelIndexArray[idx]
	print "sum(indInPixel):", sum(indInPixel)
	offsetTmp = np.zeros(cN)
	for idx1 in range(0, cN):
	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	    tmp1 = ra1[ind1][uniqueIdx1]	
	    #print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	    ind0 = obj_id0 == objInRadius[idx1]
	    if (np.sum(ind0)<2): continue
	    xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	    tmp0 = ra0[ind0][uniqueIdx0]

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	    tmp2 = ra2[ind2][uniqueIdx2]
	    #print xtmp2, xtmp1, tmp2, tmp1


	    tmp = np.concatenate((tmp0, tmp1, tmp2))
	    #tmp = tmp0
	    ytmp = (tmp-np.min(tmp))*1000*3600
	    yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2]))*1000*3600 
	    #yErrtmp = raErr0[ind0][uniqueIdx0]*1000*3600
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))  

	    if((len(xtmp)>2)):
	        num0 = len(raErr0[ind0][uniqueIdx0]) # Here num0 means the num of PS1	
	        print "num0:", num0	
	        xpm = np.zeros([4, num0])
	        for idxX in range(0, num0): # X means cross-validation
		    yErrtmp[idxX] = yErrtmp[idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		    pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		    xpm[0, idxX] = pX[0]*365.0
		    xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		    xpm[2, idxX] = (ytmp[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		    xpm[3, idxX] = pX[1]
		    yErrtmp[idxX] = yErrtmp[idxX] - 9999.0

		#idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
	        #muX = np.median(np.delete(xpm[0], np.argmax(abs(xpm[2])))) # (Here should be wrong, since removing the most import value)
		sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
		idxO = abs(xpm[2] - np.median(xpm[2]))>2.0*sigTmp
		if(sum(idxO)): 
		    idxO = np.argmax(abs(xpm[2]))
		    gmuX = xpm[0][idxO]
	            gslp = xpm[3][idxO]
		    print idx1
	        else:
		    p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		    gmuX = p[0]*365.0
		    gslp = p[1]
		    #gmuX = np.median(xpm[0]) # This case does not work very well
		    #gslp = (np.median(ytmp)-gmuX/365.0*np.median(xtmp))
	        #yErrtmp[idxO] = yErrtmp[idxO] + 9999.0
	        #p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	        #yErrtmp[idxO] = yErrtmp[idxO] - 9999.0
	        gind = gobj_id == objInRadius[idx1]
		#print objInRadius[idx1], gobj_id[gind], gra[gind], np.min(tmp), gmjd[gind], muX/365.0, slp
	        offsetTmp[idx1] = (gra[gind]-np.min(tmp))*1000*3600 - (gmjd[gind] *gmuX/365.0 + gslp)

	roffsetArray[indInPixel] = np.median(offsetTmp)  # Finally we should compare the roffsetArray and offsetArray


    ##############combining the corrected GAIA data point to fitting the PM############
    uniqueStar = np.unique(tobj_id)
    #N = 10#len(uniqueStar)
    #print "len(uniqueStar)", N
    outliers = np.array([4796157681289161166, 4796104904731190081, 4796122496917242833, 4814154487612597084, \
	4796157681289159647, 4796157681289376695]) #4796157681289159480,
    uniqueStar = outliers
    N = 7
    R2 = []

    for idx in range(0, N):
	ind1 = obj_id1 == uniqueStar[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == uniqueStar[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]

	ind2 = obj_id2 == uniqueStar[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]


	gind = gobj_id == uniqueStar[idx]
	gxtmp, guniqueIdx = np.unique(gmjd[gind], return_index=True)
	gtmp = gra[gind][guniqueIdx] - roffsetArray[gind]/(1000*3600.0) # removing the offset for each star
	print xtmp2, xtmp1, gxtmp, tmp2, tmp1, gtmp

	tmp = np.concatenate((tmp0, gtmp))#tmp1, tmp2,
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], np.array([3.0/(1000*3600.0)])))*1000*3600 #raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2], 
	xtmp = np.concatenate((xtmp0, gxtmp))#, xtmp1, xtmp2

	ind4Pal5 = tobj_id==uniqueStar[idx]
	if np.sum(ind4Pal5)>0: 
	    muPal5tmp = tmuRA[ind4Pal5][0]
	    muErrPal5tmp = tmuRAerr[ind4Pal5][0]
	else:
	    muPal5tmp = 999
	    muErrPal5tmp = 999

	if((len(xtmp)>2)):
	    num0 = len(raErr0[ind0][uniqueIdx0]) # Here num0 means the num of PS1
	    objIDtmp = uniqueStar[idx]
	    rMagtmp = mr0[ind0][0]	
	    print "num0:", num0	
	    xpm = np.zeros([4, num0])
	    for idxX in range(0, num0): # X means cross-validation
		yErrtmp[idxX] = yErrtmp[idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = (ytmp[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		yErrtmp[idxX] = yErrtmp[idxX] - 9999.0
	    
	    sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
	    idxO = abs(xpm[2] - np.median(xpm[2]))>2.0*sigTmp
	    print "sum(idxO):", sum(idxO)	    
	    if(sum(idxO)):
		idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
		gmuX = xpm[0][idxO]
	        gslp = xpm[3][idxO]
	    else:
		gmuX = np.median(xpm[0])
		gslp = (np.median(ytmp[0:num0])-gmuX/365.0*np.median(xtmp[0:num0]))

	    print "gmuX:", gmuX, xpm[0], xpm[2], sigTmp
	    #yErrtmp[idxO] = yErrtmp[idxO] + 9999.0
	    gp, gcov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    gmu = gp[0]*365 
	    gmuErr = np.sqrt(np.diag(gcov))[0]*365
	    gchi2 = np.sum((gp[0]*xtmp + gp[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
	    #yErrtmp[idxO] = yErrtmp[idxO] - 9999.0

	    p, cov = curve_fit(func, xtmp[0:-1], ytmp[0:-1], sigma=yErrtmp[0:-1], absolute_sigma=True) # remove the last point(GAIA's)
	    mu = p[0]*365 
	    muErr = np.sqrt(np.diag(cov))[0]*365
	    chi2 = np.sum((p[0]*xtmp[0:-1] + p[1] - ytmp[0:-1])**2/yErrtmp[0:-1]**2)#/(xtmp.size-1)


	    R2.append([uniqueStar[idx], rMagtmp, gchi2, gmu, gmuErr, gmuX, chi2, mu, muErr, muPal5tmp, muErrPal5tmp, roffsetArray[gind], offsetArray[gind]])


	    if(pflag):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		ax.plot(xtmp, xtmp*(muPal5tmp/365.0)+(np.median(ytmp)-muPal5tmp/365.0*np.median(xtmp)), 'g--')
		#ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')
			
		ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% gchi2 +'$(%0.2f)$'% chi2, fontsize=10)
		ax.plot(xtmp, xtmp*gmuX/365.0+gslp, 'r-')
		ax.plot(xtmp, xtmp*p[0]+p[1], 'r--')
		ax.errorbar(gmjd[gind], (gtmp-np.min(tmp))*1000*3600, 3, fmt='.', color='y')
		ax.errorbar(mjd1[ind1], (tmp1-np.min(tmp))*1000*3600, raErr1[ind1][uniqueIdx1]*1000*3600, fmt='.', color='k')
		ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		#ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_ylabel('$\delta-min(\delta)/mas$')
		ax.set_xlabel('$MJD$')
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" %gmuErr + \
		    "$(%0.2f$" %(mu)+"$\pm%0.2f)$" %muErr + \
		    ", $\Delta = %0.2f$" %(offsetArray[gind])+"$(%0.2f)$" %roffsetArray[gind] \
		    , fontsize=10) #", $\mu_{\delta}(Pal5) = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp
		ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
		figname = "fig_Pal5_GAIA_" + "%d" %idx
		plt.savefig(root+"figs/GAIA/%s_outlier.png" %figname)

    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/GAIA/R2_PS1_Pal5_GAIA.npy", R2)



##############fitting PM in different objID for Real GAIA data################
def fitPMobjGAIA2(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, raErr2, mjd2, \
    gobj_id, gra, gdec, gdecErr, gmjd, tobj_id, tmuRA, tmuRAerr, tra, tdec, pflag, root):


    ################# recovering GAIA's offset #############
    nside = 2**10 # pixel area: 0.0032784908016061202 deg^2
    phiForObj   = (gra*np.pi)/180 
    thetaForObj = (90 - gdec)* (np.pi/180) 
    pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 	    
    pixelIndexArray = np.unique(pixelIndexForObj)
    Npixel = pixelIndexArray.size
    print "Going to process %d pixels." % Npixel
    theta, phi = hp.pix2ang(nside, pixelIndexArray)
    pixelRa  = 180*phi/np.pi  # ra and dec of the centers of cells
    pixelDec = 90 - theta*180/np.pi
    #cN = 100

    roffsetArray = np.zeros(len(gdec))  # recovered offset
    roffsetErrs = np.zeros(len(gdec))
    for idx in range(0, Npixel):
        distTmp = sphdist(pixelRa[idx], pixelDec[idx], gra, gdec)
        #angSepMask = np.argsort(distTmp)

        angSepMask = (distTmp) <= (60.0/60.0)
        if(sum(angSepMask)>600):
            angSepMask = np.argsort(distTmp)
            angSepMask = angSepMask[0:600]
        elif(sum(angSepMask)<100):
            angSepMask = np.argsort(distTmp)
            angSepMask = angSepMask[0:100]

        #angSepMask = angSepMask[0:cN]
	cN = len(angSepMask)
	objInRadius = gobj_id[angSepMask]
	indInPixel = pixelIndexForObj == pixelIndexArray[idx]
	print "sum(indInPixel):", sum(indInPixel)
	offsetTmp = np.zeros(cN)
	for idx1 in range(0, cN):
	    ind1 = obj_id1 == objInRadius[idx1]
	    xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	    tmp1 = ra1[ind1][uniqueIdx1]	
	    #print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	    ind0 = obj_id0 == objInRadius[idx1]
	    if (np.sum(ind0)<2): continue
	    xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	    tmp0 = ra0[ind0][uniqueIdx0]

	    ind2 = obj_id2 == objInRadius[idx1]
	    xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	    tmp2 = ra2[ind2][uniqueIdx2]
	    #print xtmp2, xtmp1, tmp2, tmp1


	    tmp = np.concatenate((tmp0, tmp1, tmp2))
	    #tmp = tmp0
	    ytmp = (tmp-np.min(tmp))*1000*3600
	    yErrtmp = np.concatenate((raErr0[ind0][uniqueIdx0], raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2]))*1000*3600 
	    #yErrtmp = raErr0[ind0][uniqueIdx0]*1000*3600
	    xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))  

	    if((len(xtmp)>2)):
	        num0 = len(raErr0[ind0][uniqueIdx0]) # Here num0 means the num of PS1	
	        #print "num0:", num0	
	        xpm = np.zeros([4, num0])
	        for idxX in range(0, num0): # X means cross-validation
		    yErrtmp[idxX] = yErrtmp[idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		    pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		    xpm[0, idxX] = pX[0]*365.0
		    xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		    xpm[2, idxX] = (ytmp[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		    xpm[3, idxX] = pX[1]
		    yErrtmp[idxX] = yErrtmp[idxX] - 9999.0

		#idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
	        #muX = np.median(np.delete(xpm[0], np.argmax(abs(xpm[2])))) # (Here should be wrong, since removing the most import value)
		sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
		idxO = abs(xpm[2] - np.median(xpm[2]))>3.0*sigTmp
		if(sum(idxO) & (len(xtmp)>5)): 
		    idxO = np.argmax(abs(xpm[2]))
		    gmuX = xpm[0][idxO]
	            gslp = xpm[3][idxO]
		    #print idx1
	        else:
		    p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		    gmuX = p[0]*365.0
		    gslp = p[1]
		    #gmuX = np.median(xpm[0]) # This case does not work very well
		    #gslp = (np.median(ytmp)-gmuX/365.0*np.median(xtmp))
	        #yErrtmp[idxO] = yErrtmp[idxO] + 9999.0
	        #p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	        #yErrtmp[idxO] = yErrtmp[idxO] - 9999.0
	        gind = gobj_id == objInRadius[idx1]
		#Ng = sum(gind)
		#print "sum(gind):", Ng
		#print gdec[gind], gmjd[gind]
	        offsetTmp[idx1] = (gdec[gind][0]-np.min(tmp))*1000*3600 - (gmjd[gind][0] *gmuX/365.0 + gslp)

	#print "offsetTmp:", offsetTmp	
	roffsetArray[indInPixel] = np.median(offsetTmp)  # Finally we should compare the roffsetArray and offsetArray
	roffsetErrs[indInPixel] = np.sqrt((np.pi/2)/(offsetTmp.size-1)) * (0.741*(np.percentile(offsetTmp, 75) - np.percentile(offsetTmp, 25)))
	#Nval = np.histogram(offsetTmp, 40)
	#roffsetArray[indInPixel] = Nval[1][np.argmax(Nval[0])]


    ##############combining the corrected GAIA data point to fitting the PM############
    uniqueStar = np.unique(tobj_id)
    N = len(uniqueStar)
    #print "len(uniqueStar)", N
    outliers = np.array([4796157681289161166, 4796104904731190081, 4796122496917242833, 4814154487612597084, \
	4796157681289159647, 4796157681289376695]) #4796157681289159480,
    #uniqueStar = outliers
    #N = 7
    R2 = []

    for idx in range(0, N):
	ind1 = obj_id1 == uniqueStar[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == uniqueStar[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]

	ind2 = obj_id2 == uniqueStar[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]


	gind = gobj_id == uniqueStar[idx]
	gxtmp, guniqueIdx = np.unique(gmjd[gind], return_index=True)
	print gdec[gind][guniqueIdx], len(gdec[gind][guniqueIdx]), roffsetArray[gind]
	gtmp = gdec[gind][guniqueIdx] - roffsetArray[gind]/(1000*3600.0) # removing the offset for each star
	print xtmp2, xtmp1, gxtmp, tmp2, tmp1, gtmp

	tmp = np.concatenate((tmp1, tmp2, tmp0, gtmp))#
	ytmp = (tmp-np.min(tmp))*1000*3600
	gErrTmp = np.sqrt(roffsetErrs[gind]/(1000*3600.0)**2 + gdecErr[gind][guniqueIdx]**2)
	yErrtmp = np.concatenate((raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2], raErr0[ind0][uniqueIdx0], gErrTmp))*1000*3600
	xtmp = np.concatenate((xtmp1, xtmp2, xtmp0, gxtmp))#

	ind4Pal5 = tobj_id==uniqueStar[idx]
	if np.sum(ind4Pal5)>0: 
	    muPal5tmp = tmuRA[ind4Pal5][0]
	    muErrPal5tmp = tmuRAerr[ind4Pal5][0]
	else:
	    muPal5tmp = 999
	    muErrPal5tmp = 999

	if((len(xtmp)>2)):
	    s0 = len(xtmp1) + len(xtmp2)
	    num0 = len(raErr0[ind0][uniqueIdx0]) + 1 # Here num0 means the num of PS1 (+1 == including one GAIA point)
	    objIDtmp = uniqueStar[idx]
	    rMagtmp = mr0[ind0][0]	
	    print "num0:", num0	
	    xpm = np.zeros([4, num0])
	    for idxX in range(0, num0): # X means cross-validation
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] + 9999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = (ytmp[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] - 9999.0
	    
	    sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
	    idxO = abs(xpm[2] - np.median(xpm[2]))>3.0*sigTmp
	    print "sum(idxO):", sum(idxO)	    
	    if(sum(idxO) & (len(xtmp)>5)):
		idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
		gmuX = xpm[0][idxO]
	        #gslp = xpm[3][idxO]
	    else:
		gmuX = np.median(xpm[0])
		#gslp = (np.median(ytmp[0:num0])-gmuX/365.0*np.median(xtmp[0:num0]))



	    #print "gmuX:", gmuX, xpm[0], xpm[2], sigTmp
	    #yErrtmp[idxO] = yErrtmp[idxO] + 9999.0
	    gp, gcov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    gmu = gp[0]*365 
	    gmuErr = np.sqrt(np.diag(gcov))[0]*365
	    gchi2 = np.sum((gp[0]*xtmp + gp[1] - ytmp)**2/yErrtmp**2)/(xtmp.size-1)
	    #yErrtmp[idxO] = yErrtmp[idxO] - 9999.0

	    p, cov = curve_fit(func, xtmp[0:-1], ytmp[0:-1], sigma=yErrtmp[0:-1], absolute_sigma=True) # remove the last point(GAIA's)
	    mu = p[0]*365 
	    muErr = np.sqrt(np.diag(cov))[0]*365
	    chi2 = np.sum((p[0]*xtmp[0:-1] + p[1] - ytmp[0:-1])**2/yErrtmp[0:-1]**2)/(xtmp.size-1)


	    R2.append([uniqueStar[idx], rMagtmp, gchi2, gmu, gmuErr, gmuX, chi2, mu, muErr, muPal5tmp, muErrPal5tmp, roffsetArray[gind]])


	    if(pflag):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		ax.plot(xtmp, xtmp*(muPal5tmp/365.0)+(np.median(ytmp)-muPal5tmp/365.0*np.median(xtmp)), 'g--')
		#ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')
			
		ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% gchi2 +'$(%0.2f)$'% chi2, fontsize=10)
		#ax.plot(xtmp, xtmp*gmuX/365.0+gslp, 'r-')  
		ax.plot(xtmp, xtmp*(gmuX/365.25)+(np.median(ytmp)-gmuX/365.25*np.median(xtmp)), 'r-')
		ax.plot(xtmp, xtmp*p[0]+p[1], 'r--')
		ax.errorbar(gmjd[gind], (gtmp-np.min(tmp))*1000*3600, gdecErr[gind][guniqueIdx], fmt='.', color='y')
		ax.errorbar(mjd1[ind1], (tmp1-np.min(tmp))*1000*3600, raErr1[ind1][uniqueIdx1]*1000*3600, fmt='.', color='k')
		ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		#ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_ylabel('$\delta-min(\delta)/mas$')
		ax.set_xlabel('$MJD$')
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 

		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" %gmuErr + \
		    "$(%0.2f$" %(mu)+"$\pm%0.2f)$" %muErr + \
		    ", $\Delta = %0.2f$" %roffsetArray[gind] + \
		    ", $\mu_{\delta}(Pal5) = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp \
		    , fontsize=10) #
		ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
		figname = "fig_Pal5_GAIA_" + "%d" %idx
		plt.savefig(root+"figs/%s_real5_X.png" %figname)

    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/GAIA/R2_PS1_Pal5_rGAIA.npy", R2)


##############fitting PM in different objID for Real GAIA data################
def fitPMobjGAIA3(obj_id0, ra0, raErr0, mjd0, mr0, obj_id1, ra1, raErr1, mjd1, obj_id2, ra2, raErr2, mjd2, \
    gobj_id, gra, graErr, gmjd, tobj_id, tmuRA, tmuRAerr, tra, tdec, pflag, root):


    ##############combining the corrected GAIA data point to fitting the PM############
    uniqueStar = np.unique(tobj_id)
    #N = 10#len(uniqueStar)
    #print "len(uniqueStar)", N
    #outliers = np.array([4796157681289161166, 4796104904731190081, 4796122496917242833, 4814154487612597084, \
    #	4796157681289159647, 4796157681289376695]) #4796157681289159480,
    #uniqueStar = np.array([4796104904731187909])
    N = len(uniqueStar)
    R2 = []

    for idx in range(0, N):
	ind1 = obj_id1 == uniqueStar[idx]
	xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	tmp1 = ra1[ind1][uniqueIdx1]	
	#print "obj_id of SDSS:", obj_id1[ind1][0], tmp1

	ind0 = obj_id0 == uniqueStar[idx]
	if (np.sum(ind0)<2): continue
	xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	tmp0 = ra0[ind0][uniqueIdx0]

	ind2 = obj_id2 == uniqueStar[idx]
	xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	tmp2 = ra2[ind2][uniqueIdx2]


	gind = gobj_id == uniqueStar[idx]
	if (np.sum(gind)<1): continue
	gxtmp, guniqueIdx = np.unique(gmjd[gind], return_index=True)
	gtmp = gra[gind][guniqueIdx]
	#print xtmp2, xtmp1, gxtmp, tmp2, tmp1, gtmp

	tmp = np.concatenate((tmp1,tmp2, tmp0, gtmp))# 
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = np.concatenate((raErr1[ind1][uniqueIdx1], raErr2[ind2][uniqueIdx2], raErr0[ind0][uniqueIdx0], graErr[gind][guniqueIdx]))*1000*3600 # 
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))#  

	ind4Pal5 = tobj_id==uniqueStar[idx]
	if np.sum(ind4Pal5)>0: 
	    muPal5tmp = tmuRA[ind4Pal5][0]
	    muErrPal5tmp = tmuRAerr[ind4Pal5][0]
	else:
	    muPal5tmp = 999
	    muErrPal5tmp = 999

	if((len(xtmp)>2)):
	    s0 = len(xtmp1) + len(xtmp2) # 
	    num0 = len(raErr0[ind0][uniqueIdx0]) + len(gtmp) # Here num0 means the num of PS1
	    objIDtmp = uniqueStar[idx]
	    rMagtmp = mr0[ind0][0]	
	    print "num0:", num0	
	    xpm = np.zeros([5, num0])
	    for idxX in range(0, num0): # X means cross-validation
		#print  yErrtmp[s0+idxX]
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] + 99999.0 #loop ONLY for PS1 data, so "num0 + idxX"
		pX, covX = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		xpm[0, idxX] = pX[0]*365.0
		xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.0
		xpm[2, idxX] = abs(ytmp[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
		xpm[3, idxX] = pX[1]
		xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - np.delete(ytmp, s0+idxX))**2/np.delete(yErrtmp, s0+idxX)**2)/(xtmp.size-2)
		yErrtmp[s0+idxX] = yErrtmp[s0+idxX] - 99999.0
	    
	    #sigTmp = (0.741*(np.percentile(xpm[2], 75) - np.percentile(xpm[2], 25)))
	    #idxO = abs(xpm[2] - np.median(xpm[2]))>2.0*sigTmp
	    #print "sum(idxO):", sum(idxO), xpm[4], xpm[0]	    
	    '''if(sum(idxO) & (len(xtmp)>5)):
		idxO = np.argmax(abs(xpm[2]))  # the max delta are probably a outlier
		gmuX = xpm[0][idxO]
	    else:
		gmuX = np.median(xpm[0])'''
	    
	    gp, gcov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
	    gmu = gp[0]*365 
	    gmuErr = np.sqrt(np.diag(gcov))[0]*365
	    gchi2 = np.sum((gp[0]*xtmp + gp[1] - ytmp)**2/yErrtmp**2)/(xtmp.size-1)
	    #print len(xpm[4]), xpm[4], gchi2
	    chi2Tmp = np.concatenate((xpm[4], np.array([gchi2])))
	    gmuTmp = np.concatenate((xpm[0], np.array([gmu])))
	    gslpTmp = np.concatenate((xpm[3], np.array([gp[1]])))
 	    idxO = np.argmin(chi2Tmp) 
	    gmuX = gmuTmp[idxO]
	    gslp = gslpTmp[idxO]
	    gchi2 = np.min(chi2Tmp)

	    #p, cov = curve_fit(func, xtmp[0:-1], ytmp[0:-1], sigma=yErrtmp[0:-1], absolute_sigma=True) # remove the last point(GAIA's)
	    #mu = p[0]*365 
	    #muErr = np.sqrt(np.diag(cov))[0]*365
	    #chi2 = np.sum((p[0]*xtmp[0:-1] + p[1] - ytmp[0:-1])**2/yErrtmp[0:-1]**2)/(xtmp.size-2)

	    mu = xpm[0][num0-1]
	    muErr = xpm[1][num0-1]
	    chi2 = xpm[4][num0-1]

	    R2.append([uniqueStar[idx], rMagtmp, gchi2, gmu, gmuErr, gmuX, chi2, mu, muErr, muPal5tmp, muErrPal5tmp])


	    if(pflag):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		ax.plot(xtmp, xtmp*(muPal5tmp/365.0)+(np.median(ytmp)-muPal5tmp/365.0*np.median(xtmp)), 'g--')
		#ax.plot(xtmp, xtmp*(muPS1tmp/365.0)+(np.median(ytmp)-muPS1tmp/365.0*np.median(xtmp)), 'k--')
			
		ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% gchi2 +'$(%0.2f)$'% chi2, fontsize=10)
		ax.plot(xtmp, xtmp*gmuX/365.0+gslp, 'r-')  
		#ax.plot(xtmp, xtmp*(gmuX/365.25)+(np.median(ytmp)-gmuX/365.25*np.median(xtmp)), 'r-')
		#ax.plot(xtmp, xtmp*p[0]+p[1], 'r--')
		ax.errorbar(gmjd[gind], (gtmp-np.min(tmp))*1000*3600, graErr[gind][guniqueIdx], fmt='.', color='y')
		#ax.errorbar(mjd1[ind1], (tmp1-np.min(tmp))*1000*3600, raErr1[ind1][uniqueIdx1]*1000*3600, fmt='.', color='k')
		ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		#ax.set_ylabel('$\eta-min(\eta)/mas$')
		ax.set_ylabel('$\delta-min(\delta)/mas$')
		ax.set_xlabel('$MJD$')
		ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100)

		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    %np.round(rMagtmp,decimals=2), fontsize=10)
		ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu_{\delta} = %0.2f$" %(gmuX)+"$\pm%0.2f$" %gmuErr + \
		    "$(%0.2f$" %(mu)+"$\pm%0.2f)$" %muErr + \
		    ", $\mu_{\delta}(Pal5) = %0.2f$" %(muPal5tmp)+"$\pm%0.2f$" %muErrPal5tmp, \
		    fontsize=10)
		ax.set_title('Comparision of Proper Motions Using Different Methods', fontsize=12)
		figname = "fig_Pal5_GAIA_" + "%d" %idx
		plt.savefig(root+"figs/%s_real_parallel_GP1_100_tst6.png" %figname)

    if(not pflag):
    	R2 = np.array(R2,dtype=object)
    	np.save(root+"figs/R2_Pal5_rGAIA_GP1_80.npy", R2)

##############fitting PM in different objID################
#fittingPM is same as fittingPM2, but fittingPM is too slow, because operating directly on table
def fittingPM(table0, table1, table2, magMask, pflag, chunkNo, ra0, dec0, root): 
    uniqueStar = np.unique(table0['obj_id'])
    tnum = uniqueStar.size
    pm = []
    if(pflag):
	outliers = np.array([-4997992831776515880, -4799869632544364352, -4763770466781189284, -4763770466781242308, \
	    -4799869632544364323])#
	objID_curr = uniqueStar#outliers
	N = 10#len(uniqueStar)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	flag = 0  #record the SDSS 
	objIDtmp = objID_curr[idx]
	ind0 = table0['obj_id'] == objIDtmp
	if (np.sum(ind0)<2): continue
	#xtmp0, uniqueIdx0 = np.unique(table0['mjd'][ind0], return_index=True)
	xtmp0 = table0['mjd'][ind0]
	xitmp0, etatmp0, status = s2t.ds2tp(table0['ra'][ind0], table0['dec'][ind0], ra0, dec0)

	ind1 = table1.get_where_list('(obj_id == objIDtmp)')
	#xtmp1, uniqueIdx1 = np.unique(table1.col('mjd')[ind1], return_index=True)
	xtmp1 = table1.col('mjd')[ind1]
	xitmp1, etatmp1, status = s2t.ds2tp(table1.col('ra')[ind1], table1.col('dec')[ind1], ra0, dec0)
	if(ind1.size>0): flag = flag + 10	

	ind2 = table2.get_where_list('(obj_id == objIDtmp)')
	#xtmp2, uniqueIdx2 = np.unique(table2.col('mjd')[ind2], return_index=True)
	xtmp2 = table2.col('mjd')[ind2]
	xitmp2, etatmp2, status = s2t.ds2tp(table2.col('ra')[ind2], table2.col('dec')[ind2], ra0, dec0)
	if(ind2.size>0): flag = flag + 5

	xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))# 
	etatmp = np.concatenate((etatmp0, etatmp1, etatmp2))# 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((table0['raErr'][ind0], table1.col('raErr')[ind1], \
	    table2.col('raErr')[ind2]))*1000*3600 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((table0['decErr'][ind0], table1.col('decErr')[ind1], \
	    table2.col('decErr')[ind2]))*1000*3600

	xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 

	if((len(xtmp)>2)):
		rMagtmp = table0['mr'][ind0][0]
		num0 = np.sum(ind0)	
		print "starIdx, tnum, num0, chunkNo:", idx, tnum, num0,	chunkNo
		xpm = np.zeros([6, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[idxX] = yErrtmpXI[idxX] + 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] + 9999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
			xpm[3, idxX] = pXeta[0]*365.25
			xpm[4, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[5, idxX] = (ytmpETA[idxX] - (pXeta[0]*xtmp[idxX] + pXeta[1]))
			yErrtmpXI[idxX] = yErrtmpXI[idxX] - 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] - 9999.0
		muXxi = np.median(xpm[0])
		muErrXxi = (0.741*(np.percentile(xpm[0], 75) - np.percentile(xpm[0], 25)))
		muXeta = np.median(xpm[3])
		muErrXeta = (0.741*(np.percentile(xpm[3], 75) - np.percentile(xpm[3], 25)))


		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		variance = np.var(ytmpXI)
		residuals = np.var(ytmpXI - pXI[0]*xtmp + pXI[1])
		Rtmp = residuals/variance		
		R2XI = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)#/(xtmp.size-1)

		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		variance = np.var(ytmpETA)
		residuals = np.var(ytmpETA - pETA[0]*xtmp + pETA[1])
		Rtmp = residuals/variance		
		R2ETA = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)

		#print objIDtmp, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, chi2ETA, muETA, muErrETA, muXeta, muErrXeta
		pm.append([objIDtmp, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag])
		if(pflag):
			p = pETA
			mu = muETA
			muErr = muErrETA
			yErrtmp = yErrtmpETA
			ytmp = ytmpETA
			muX = muXeta
			muErrX = muErrXeta
			R2 = R2ETA
			chi2 = chi2ETA

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_dec_" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()
    if(not pflag):
    	pm = np.array(pm,dtype=object)
    	np.save(root+"figs/PMs"+str(chunkNo)+".npy", pm)

##############fitting PM in different objID################
def fittingPM2(table0, table1, table2, magMask, pflag, chunkNo, root):

    obj_id0 = table0['obj_id']
    ra0 = table0['ra']
    raErr0 = table0['raErr']
    dec0 = table0['dec']
    decErr0 = table0['decErr']
    mjd0 = table0['mjd']
    mr0 = table0['mr']

    CRA = np.median(ra0)
    CDEC = np.median(dec0)

    if(table1==0):
        obj_id1 = np.zeros(10) - 999
        ra1 = np.zeros(10)
        raErr1 = np.zeros(10)
        dec1 = np.zeros(10)
        decErr1 = np.zeros(10)
        mjd1 = np.zeros(10)
	print "No coverage in SDSS this chunk"	
    else:
        obj_id1 = table1.col('obj_id')
        ra1 = table1.col('ra')
        raErr1 = table1.col('raErr')
        dec1 = table1.col('dec')
        decErr1 = table1.col('decErr')
        mjd1 = table1.col('mjd')

    obj_id2 = table2.col('obj_id')
    ra2 = table2.col('ra')
    raErr2 = table2.col('raErr')
    dec2 = table2.col('dec')
    decErr2 = table2.col('decErr')
    mjd2 = table2.col('mjd')

    uniqueStar = np.unique(obj_id0)
    tnum = uniqueStar.size
    pm = []
    if(pflag):
	outliers = np.array([4796104904731190081, -4997992831776515880, -4799869632544364352, -4763770466781189284, -4763770466781242308, \
	    -4799869632544364323])#
	objID_curr = outliers#uniqueStar#
	N = 1#len(uniqueStar)
    else:
	N=len(uniqueStar)
	objID_curr = uniqueStar
    for idx in range(0, N):
	flag = 0  #record the SDSS 
	objIDtmp = objID_curr[idx]
	ind0 = obj_id0 == objIDtmp
	if (np.sum(ind0)<2): continue
	#xtmp0, uniqueIdx0 = np.unique(mjd0[ind0], return_index=True)
	xtmp0 = mjd0[ind0]
	xitmp0, etatmp0, status = s2t.ds2tp(ra0[ind0], dec0[ind0], CRA, CDEC)

	ind1 = obj_id1 == objIDtmp
	#xtmp1, uniqueIdx1 = np.unique(mjd1[ind1], return_index=True)
	xtmp1 = mjd1[ind1]
	xitmp1, etatmp1, status = s2t.ds2tp(ra1[ind1], dec1[ind1], CRA, CDEC)
	if(np.sum(ind1)>0): flag = flag + 10	

	ind2 = obj_id2 == objIDtmp
	#xtmp2, uniqueIdx2 = np.unique(mjd2[ind2], return_index=True)
	xtmp2 = mjd2[ind2]
	xitmp2, etatmp2, status = s2t.ds2tp(ra2[ind2], dec2[ind2], CRA, CDEC)
	if(np.sum(ind2)>0): flag = flag + 5

	xitmp = np.concatenate((xitmp0, xitmp1, xitmp2))# 
	etatmp = np.concatenate((etatmp0, etatmp1, etatmp2))# 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr0[ind0], raErr1[ind1], \
	    raErr2[ind2]))*1000*3600 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr0[ind0], decErr1[ind1], \
	    decErr2[ind2]))*1000*3600

	xtmp = np.concatenate((xtmp0, xtmp1, xtmp2))# 


	if((len(xtmp)>2)):
		rMagtmp = mr0[ind0][0]
		num0 = np.sum(ind0)	
		print "starIdx, tnum, num0, chunkNo:", idx, tnum, num0,	chunkNo
		xpm = np.zeros([6, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[idxX] = yErrtmpXI[idxX] + 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] + 9999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
			xpm[3, idxX] = pXeta[0]*365.25
			xpm[4, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[5, idxX] = (ytmpETA[idxX] - (pXeta[0]*xtmp[idxX] + pXeta[1]))
			yErrtmpXI[idxX] = yErrtmpXI[idxX] - 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] - 9999.0
		muXxi = np.median(xpm[0])
		muErrXxi = (0.741*(np.percentile(xpm[0], 75) - np.percentile(xpm[0], 25)))
		muXeta = np.median(xpm[3])
		muErrXeta = (0.741*(np.percentile(xpm[3], 75) - np.percentile(xpm[3], 25)))


		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		variance = np.var(ytmpXI)
		residuals = np.var(ytmpXI - pXI[0]*xtmp + pXI[1])
		Rtmp = residuals/variance		
		R2XI = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)#/(xtmp.size-1)

		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		variance = np.var(ytmpETA)
		residuals = np.var(ytmpETA - pETA[0]*xtmp + pETA[1])
		Rtmp = residuals/variance		
		R2ETA = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)

		#print objIDtmp, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, chi2ETA, muETA, muErrETA, muXeta, muErrXeta
		pm.append([objIDtmp, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, \
		    chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag, ra0[ind0][0], dec0[ind0][0]])
		if(pflag):
			p = pETA
			mu = muETA
			muErr = muErrETA
			yErrtmp = yErrtmpETA
			ytmp = ytmpETA
			muX = muXeta
			muErrX = muErrXeta
			R2 = R2ETA
			chi2 = chi2ETA

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %objIDtmp, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_dec_" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()
    if(not pflag):
    	pm = np.array(pm,dtype=object)
    	np.save(root+"figs/PM_radec_"+str(chunkNo)+".npy", pm)

##############fitting PM for each star parallelly################
'''def fittingPM3(packParameterList):
    idx, uniqueID, obj_id0, ra0, dec0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, pflag, chunkNo, CRA, CDEC, root = packParameterList

    flag = 0  #record the PS1 
    ind0 = obj_id0 == uniqueID
    #print 'ind0', ind0
    if(np.sum(ind0)>2):
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
	    #raErr2[ind2]))*1000*3600 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	#yErrtmpETA = 0.000000001+np.concatenate((decErr1[ind1], \
	    #decErr2[ind2]))*1000*3600

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
    yErrtmpXI = np.zeros([np.sum(ind0)])
    yErrtmpETA = np.zeros([np.sum(ind0)])
    print 'len', len(xtmp)
    if((len(xtmp)>0)):
		rMagtmp = mr0[ind0][0]
		print(1)
		num0 = np.sum(ind0)	
		#print "starIdx, num0, chunkNo:", idx, num0,chunkNo
		xpm = np.zeros([6, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[idxX] = yErrtmpXI[idxX] + 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] + 9999.0
			#print yErrtmpXI
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365
			xpm[2, idxX] = (ytmpXI[idxX] - (pX[0]*xtmp[idxX] + pX[1]))
			xpm[3, idxX] = pXeta[0]*365
			xpm[4, idxX] = np.sqrt(np.diag(covXeta))[0]*365
			xpm[5, idxX] = (ytmpETA[idxX] - (pXeta[0]*xtmp[idxX] + pXeta[1]))
			yErrtmpXI[idxX] = yErrtmpXI[idxX] - 9999.0
			yErrtmpETA[idxX] = yErrtmpETA[idxX] - 9999.0
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
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			#ax.plot(xtmp, xtmp*(muX/365)+(np.median(ytmp)-muX/365*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'b-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta$(mas)')
			ax.set_xlabel('MJD')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			#ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    #%muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  #$flag = %d$" %(flag), fontsize=10)
			#ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    #%muErr + ", $flag = %d$" %(flag), fontsize=10)

			ax.set_title('Proper Motions with Only the Simulated CSST Survey', fontsize=12)
			figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s_Mock.png" %figname)
			#plt.show()
    return uniqueID, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, \
	    chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag, ra0[ind0][0], dec0[ind0][0]
    #else:
	#return idx, uniqueID'''
    


def fittingPM3(packParameterList):
    idx, uniqueID, obj_id0, ra0, dec0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, pflag, chunkNo, CRA, CDEC, root = packParameterList

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
    yErrtmpXI = np.ones([np.sum(ind0)])
    yErrtmpETA = np.ones([np.sum(ind0)])
    errbar = np.zeros([np.sum(ind0)])
    #print 'len', len(xtmp)
    if((len(xtmp)>0)):
        rMagtmp = mr0[ind0][0]

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

            #print xtmp, ytmpXI, yErrtmpXI
            pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
            muXI = pXI[0]*365
            #print muXI
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
                ax.errorbar(xtmp, ytmp, errbar, fmt='.', color='red')
			#ax.plot(xtmp, xtmp*(muX/365)+(np.median(ytmp)-muX/365*np.median(xtmp)), 'b--')
                #ax.text(np.min(xtmp)+50, np.max(ytmp)*0.72,'PM = %0.2f'% ytmp[0], fontsize=10)
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
                        %muErr +", flag = %d" %(flag), fontsize=10)
                #print(muErr)
                ax.set_title('Proper Motions with M31', fontsize=12)
                figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
                plt.savefig(root+"figs/%s_Mock.png" %figname)
			#plt.show()
			#plt.show()
    return uniqueID, rMagtmp, R2XI, chi2XI, muXI, muErrXI, muXxi, muErrXxi, R2ETA, \
    	    chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag, ra0[ind0][0], dec0[ind0][0]
    #else:
	#return idx, uniqueID'''


##############fitting PM for each star parallelly, + GAIA################
def fittingPM4(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	pflag, chunkNo, CRA, CDEC, root = packParameterList

    flag = 0  #record the PS1 
    ind0 = obj_id0 == uniqueID
    if(np.sum(ind0)>2):

	xtmp0 = mjd0[ind0]
	CRA = np.median(ra0[ind0])
 	CDEC = np.median(dec0[ind0])
	xitmp0, etatmp0, status = s2t.ds2tp(ra0[ind0], dec0[ind0], CRA, CDEC)

	ind1 = obj_id1 == uniqueID
	xtmp1 = mjd1[ind1]
	xitmp1, etatmp1, status = s2t.ds2tp(ra1[ind1], dec1[ind1], CRA, CDEC)
	if(np.sum(ind1)>0): flag = flag + 10	

	ind2 = obj_id2 == uniqueID
	xtmp2 = mjd2[ind2]
	xitmp2, etatmp2, status = s2t.ds2tp(ra2[ind2], dec2[ind2], CRA, CDEC)
	if(np.sum(ind2)>0): flag = flag + 5

	gind = gobj_id == uniqueID
	gxtmp = gmjd[gind]
	gxitmp, getatmp, status = s2t.ds2tp(gra[gind], gdec[gind], CRA, CDEC)
	if(np.sum(gind)>0): flag = flag + 20

	xitmp = np.concatenate((xitmp1, xitmp2, xitmp0, gxitmp))
	etatmp = np.concatenate((etatmp1, etatmp2, etatmp0, getatmp))

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1[ind1], \
	    raErr2[ind2], raErr0[ind0], graErr[gind]))*1000*3600 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1[ind1], \
	    decErr2[ind2], decErr0[ind0], gdecErr[gind]))*1000*3600
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))#  
	if((len(xtmp)>2)):
		rMagtmp = mr0[ind0][0]
		num0 = np.sum(ind0)  + len(gxtmp) #X-validate PS1 + GAIA points
		s0 = len(xtmp1) + len(xtmp2) #	
		print "starIdx, num0, chunkNo:", idx, num0,chunkNo
		xpm = np.zeros([10, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] + 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] + 999999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
			xpm[3, idxX] = pX[1]
			xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - \
			    np.delete(ytmpXI, s0+idxX))**2/np.delete(yErrtmpXI, s0+idxX)**2)/(xtmp.size-2)

			xpm[5, idxX] = pXeta[0]*365.25
			xpm[6, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[7, idxX] = (ytmpETA[s0+idxX] - (pXeta[0]*xtmp[s0+idxX] + pXeta[1]))
			xpm[8, idxX] = pXeta[1]
			xpm[9, idxX] = np.sum((pXeta[0]*np.delete(xtmp, s0+idxX) + pXeta[1] - \
			    np.delete(ytmpETA, s0+idxX))**2/np.delete(yErrtmpETA, s0+idxX)**2)/(xtmp.size-2)

			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] - 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] - 999999.0

		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		variance = np.var(ytmpXI)
		residuals = np.var(ytmpXI - pXI[0]*xtmp + pXI[1])
		#Rtmp = residuals/variance		
		#R2XI = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[4], np.array([chi2XI])))
	        muTmp = np.concatenate((xpm[0], np.array([muXI])))
	        slpTmp = np.concatenate((xpm[3], np.array([pXI[1]])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXxi = muTmp[idxO]
	        slpxi = slpTmp[idxO]
	        chi2XI = np.min(chi2Tmp)
		muErrXxi = (0.741*(np.percentile(xpm[0], 75) - np.percentile(xpm[0], 25))) # This value is useless


		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		variance = np.var(ytmpETA)
		residuals = np.var(ytmpETA - pETA[0]*xtmp + pETA[1])
		#Rtmp = residuals/variance		
		#R2ETA = np.round(np.abs(1-Rtmp), decimals=2) 
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[9], np.array([chi2ETA])))
	        muTmp = np.concatenate((xpm[5], np.array([muETA])))
	        slpTmp = np.concatenate((xpm[8], np.array([pETA[1]])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXeta = muTmp[idxO]
	        slpeta = slpTmp[idxO]
	        chi2ETA = np.min(chi2Tmp)
		muErrXeta = (0.741*(np.percentile(xpm[5], 75) - np.percentile(xpm[5], 25))) # This value is useless

		if(pflag):
			p = pETA#pXI
			mu = muETA#muXI
			muErr = muErrETA#muErrXI
			yErrtmp = yErrtmpETA#yErrtmpXI
			ytmp = ytmpETA#ytmpXI
			muX = muXeta#muXxi
			muErrX = muErrXeta#muErrXxi
			#R2 = R2ETA#R2XI
			chi2 = chi2ETA#chi2XI

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			#ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    #%muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    #",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s_Mock.png" %figname)
			#plt.show()
        return uniqueID, rMagtmp, chi2XI, muXI, muErrXI, muXxi, muErrXxi, \
	    chi2ETA, muETA, muErrETA, muXeta, muErrXeta, xpm, flag, ra0[ind0][0], dec0[ind0][0]
    else:
	return idx, uniqueID

##############fitting PM for each star parallelly, + GAIA, speed up based on fittingPM4 ################
def fittingPM5(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
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

	xitmp = np.concatenate((ra1, ra2, ra0, gra))#   #
	etatmp = np.concatenate((dec1, dec2, dec0, gdec))#  # 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0, graErr))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0, gdecErr))*1000*3600 # #
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))#  
	#print xtmp
	#print "preparing time:", time()-t0
	if((len(xtmp)>2)):
		rMagtmp = mr0[0]
		#num0 = len(xtmp0)  + len(gxtmp) #X-validate PS1 + GAIA points
		#s0 = len(xtmp1) + len(xtmp2) #	
		num0 = len(xtmp0)  + len(gxtmp) + len(xtmp1) #X-validate PS1 + GAIA + SDSS points
		s0 = len(xtmp2) #	
		#print "starIdx, num0, chunkNo:", idx, num0,chunkNo
		xpm = np.zeros([10, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] + 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] + 999999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
			xpm[3, idxX] = pX[1]
			xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - \
			    np.delete(ytmpXI, s0+idxX))**2/np.delete(yErrtmpXI, s0+idxX)**2)/(xtmp.size-2)

			xpm[5, idxX] = pXeta[0]*365.25
			xpm[6, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[7, idxX] = (ytmpETA[s0+idxX] - (pXeta[0]*xtmp[s0+idxX] + pXeta[1]))
			xpm[8, idxX] = pXeta[1]
			xpm[9, idxX] = np.sum((pXeta[0]*np.delete(xtmp, s0+idxX) + pXeta[1] - \
			    np.delete(ytmpETA, s0+idxX))**2/np.delete(yErrtmpETA, s0+idxX)**2)/(xtmp.size-2)

			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] - 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] - 999999.0

		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25 # fit with all the data
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[4], np.array([chi2XI])))
	        muTmp = np.concatenate((xpm[0], np.array([muXI])))
	        slpTmp = np.concatenate((xpm[3], np.array([pXI[1]])))
	        errTmp = np.concatenate((xpm[1], np.array([muErrXI])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXxi = muTmp[idxO] # it should be the best one
	        slpxi = slpTmp[idxO]
		muErrXxi = errTmp[idxO]
	        chi2XI = np.min(chi2Tmp)
		
		muXIog = xpm[0][num0-1] # the last one should be the fit without GAIA
		muErrXIog = xpm[1][num0-1]

		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25 # fit with all the data
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[9], np.array([chi2ETA])))
	        muTmp = np.concatenate((xpm[5], np.array([muETA])))
	        slpTmp = np.concatenate((xpm[8], np.array([pETA[1]])))
	        errTmp = np.concatenate((xpm[6], np.array([muErrETA])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXeta = muTmp[idxO] # it should be the best
	        slpeta = slpTmp[idxO]
		muErrXeta = errTmp[idxO]
	        chi2ETA = np.min(chi2Tmp)

		muETAog = xpm[5][num0-1] # the last one should be the fit without GAIA
		muErrETAog = xpm[6][num0-1]

		if(pflag):
			p = pETA#pXI
			mu = muETA#muXI
			muErr = muErrETA#muErrXI
			yErrtmp = yErrtmpETA#yErrtmpXI
			ytmp = ytmpETA#ytmpXI
			muX = muXeta#muXxi
			muErrX = muErrXeta#muErrXxi
			#R2 = R2ETA#R2XI
			chi2 = chi2ETA#chi2XI

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()

        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
        return uniqueID, ra, dec, rMagtmp, muXI, muErrXI, muXxi, muErrXxi, muXIog, muErrXIog, chi2XI, \
	    muETA, muErrETA, muXeta, muErrXeta, muETAog, muErrETAog, chi2ETA, len(xtmp), flag
    else:
	return idx, uniqueID


##############fitting PM to predict each GAIA star, and get the positional offset between predict and calibrated position, parallelly ################
def fittingPM6(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
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
		rMagtmp = mr0[0]


		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25 # fit with all the data
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-1) # we should use the reduced chi2



		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25 # fit with all the data
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-1) # we should use the reduced chi2


	        offsetXi = ((gra[0]-np.min(xitmp))*1000*3600 - (gxtmp[0] * pXI[0] + pXI[1]))
	        offsetEta = ((gdec[0]-np.min(etatmp))*1000*3600 - (gxtmp[0] * pETA[0] + pETA[1]))


        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(gra[0]), np.radians(gdec[0]), CRA, CDEC)
        #################################################################################
        return uniqueID, ra, dec, rMagtmp, muXI, muErrXI, offsetXi, np.min(xitmp), pXI[0], pXI[1], chi2XI, \
	    muETA, muErrETA, offsetEta, np.min(etatmp), pETA[0], pETA[1], chi2ETA, len(xtmp), flag
    else:
	return idx, uniqueID


##############fitting PM for each star parallelly, + GAIA, to invastigate the PM offset pattern based on fittingPM5(similar fittingPM4) ################
def fittingPM7(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
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

	ratmp = np.concatenate((ra1, ra2, ra0, gra))
	dectmp = np.concatenate((dec1, dec2, dec0, gdec))

	#########################(RA, DEC) ==> (xi, eta)#################################
	# use the median position of fitted star, to transform the (RA, DEC)/degree into (xi, eta)/degree
	CRA = np.median(ratmp)
 	CDEC = np.median(dectmp)
	xitmp, etatmp, status = s2t.ds2tp(ratmp, dectmp, CRA, CDEC)
	################################################################################


	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0, graErr))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0, gdecErr))*1000*3600 # #
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))#  
	#print xtmp
	#print "preparing time:", time()-t0
	if((len(xtmp)>2)):
		rMagtmp = mr0[0]
		#num0 = len(xtmp0)  + len(gxtmp) #X-validate PS1 + GAIA points
		#s0 = len(xtmp1) + len(xtmp2) #	
		num0 = len(xtmp0)  + len(gxtmp) + len(xtmp1) + len(xtmp2) #X-validate PS1 + GAIA + SDSS + 2MASS points
		s0 = 0#len(xtmp2) #	
		#print "starIdx, num0, chunkNo:", idx, num0,chunkNo
		xpm = np.zeros([10, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] + 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] + 999999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
			xpm[3, idxX] = pX[1]
			xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - \
			    np.delete(ytmpXI, s0+idxX))**2/np.delete(yErrtmpXI, s0+idxX)**2)/(xtmp.size-3)

			xpm[5, idxX] = pXeta[0]*365.25
			xpm[6, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[7, idxX] = (ytmpETA[s0+idxX] - (pXeta[0]*xtmp[s0+idxX] + pXeta[1]))
			xpm[8, idxX] = pXeta[1]
			xpm[9, idxX] = np.sum((pXeta[0]*np.delete(xtmp, s0+idxX) + pXeta[1] - \
			    np.delete(ytmpETA, s0+idxX))**2/np.delete(yErrtmpETA, s0+idxX)**2)/(xtmp.size-3)

			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] - 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] - 999999.0

		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25 # fit with all the data
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-2) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[4], np.array([chi2XI])))
	        muTmp = np.concatenate((xpm[0], np.array([muXI])))
	        slpTmp = np.concatenate((xpm[3], np.array([pXI[1]])))
	        errTmp = np.concatenate((xpm[1], np.array([muErrXI])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXxi = muTmp[idxO] # it should be the best one
	        slpxi = slpTmp[idxO]
		muErrXxi = errTmp[idxO]
	        chi2XI = np.min(chi2Tmp)
		
		if(len(gxtmp)>0):
			muXIog = xpm[0][num0-1] # the last one should be the fit without GAIA
			muErrXIog = xpm[1][num0-1]
		else:
			muXIog = muXI
			muErrXIog = muErrXI

		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25 # fit with all the data
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-2) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[9], np.array([chi2ETA])))
	        muTmp = np.concatenate((xpm[5], np.array([muETA])))
	        slpTmp = np.concatenate((xpm[8], np.array([pETA[1]])))
	        errTmp = np.concatenate((xpm[6], np.array([muErrETA])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXeta = muTmp[idxO] # it should be the best
	        slpeta = slpTmp[idxO]
		muErrXeta = errTmp[idxO]
	        chi2ETA = np.min(chi2Tmp)

		if(len(gxtmp)>0):
			muETAog = xpm[5][num0-1] # the last one should be the fit without GAIA
			muErrETAog = xpm[6][num0-1]
		else:
			muETAog = muETA
			muErrETAog = muErrETA

		####### ONLY PS1 PM#########
		sp0 = len(xtmp1) + len(xtmp2)
		onp = len(xtmp0)
		pXIps1, covXIps1 = curve_fit(func, xtmp0, ytmpXI[sp0:sp0+onp], sigma=yErrtmpXI[sp0:sp0+onp], absolute_sigma=True)
		muXIps1 = pXIps1[0]*365.25 # fit with all the data
		muErrXIps1 = np.sqrt(np.diag(covXIps1))[0]*365.25
		chi2XIps1 = np.sum((pXIps1[0]*xtmp0 + pXIps1[1] - ytmpXI[sp0:sp0+onp])**2/yErrtmpXI[sp0:sp0+onp]**2)/(onp-2) 

		pETAps1, covETAps1 = curve_fit(func, xtmp0, ytmpETA[sp0:sp0+onp], sigma=yErrtmpETA[sp0:sp0+onp], absolute_sigma=True)
		muETAps1 = pETAps1[0]*365.25 # fit with all the data
		muErrETAps1 = np.sqrt(np.diag(covETAps1))[0]*365.25
		chi2ETAps1 = np.sum((pETAps1[0]*xtmp0 + pETAps1[1] - ytmpETA[sp0:sp0+onp])**2/yErrtmpETA[sp0:sp0+onp]**2)/(onp-2) 


		if(pflag):
			p = pXI#pETA#
			mu = muXI#muETA#
			muErr = muErrXI#muErrETA#
			yErrtmp = yErrtmpETA#yErrtmpXI
			ytmp = ytmpXI#ytmpETA#
			muX = muXxi#muXeta#
			muErrX = muErrXxi#muErrXeta#
			#R2 = R2ETA#R2XI
			chi2 = chi2XI#chi2ETA#

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			#ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_XI" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()

        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        #ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
        return uniqueID, ra0[0], dec0[0], rMagtmp, muXI, muErrXI, muXxi, muErrXxi, muXIog, muErrXIog, chi2XI, \
	    muETA, muErrETA, muXeta, muErrXeta, muETAog, muErrETAog, chi2ETA, \
	    muXIps1, muErrXIps1, chi2XIps1, muETAps1, muErrETAps1, chi2ETAps1, len(xtmp0), len(xtmp), flag
    else:
	return idx, uniqueID



def lnprob(x, xtmp, ytmp, yErrtmp, ytmpETA, yErrtmpETA, gpmtmp):
    if((abs(x[0])<1.5) & (abs(x[2])<1.5)):
        p1 = -0.5*(ytmp-(x[0]*xtmp + x[1]))**2/(yErrtmp**2)-0.5*np.log(yErrtmp**2)
        p3 = -0.5*(ytmpETA-(x[2]*xtmp + x[3]))**2/(yErrtmpETA**2)-0.5*np.log(yErrtmpETA**2)
        
        if((gpmtmp[0]<9998) & (gpmtmp[1]<9998)):
            p2 = -0.5*(x[0]-(gpmtmp[0]/365.24))**2/((gpmtmp[1]/365.24)**2)-0.5*np.log((gpmtmp[1]/365.24)**2)
        else:
            p2=0
        if((gpmtmp[2]<9998) & (gpmtmp[3]<9998)):
            p4 = -0.5*(x[2]-(gpmtmp[2]/365.24))**2/((gpmtmp[3]/365.24)**2)-0.5*np.log((gpmtmp[3]/365.24)**2)
        else:
            p4=0
        lnlkh = p2 + np.sum(p1[np.where(p1>-1e10)]) + np.sum(p3[np.where(p3>-1e10)]) + p4	
        #print p2, p4, p1, p3, lnlkh
        return lnlkh
    else:
        return -np.inf

##############fitting PM for each star parallelly, + GAIA DR2, with MCMC ################
def fittingPM_mcmc(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
        obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
        dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
        gpmra, gpmdec, gpmraErr, gpmdecErr, \
        pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
    ind0 = obj_id0 == uniqueID
    xtmp0 = mjd0[ind0] 
    if(idx>10): pflag=False
    if(len(xtmp0)>2):
	ra0 = ra0[ind0]
	raErr0 = raErr0[ind0]
	dec0 = dec0[ind0]
	decErr0 = decErr0[ind0]
	mr0 = mr0[ind0]
	#print uniqueID, ra0, dec0, mr0

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
	gpmraErr = gpmraErr[gind]
	gpmdecErr = gpmdecErr[gind]
	if(len(gxtmp)>0): flag = flag + 20

	ratmp = np.concatenate((ra1, ra2, ra0, gra))
	dectmp = np.concatenate((dec1, dec2, dec0, gdec))
	#print ratmp, flag

	#########################(RA, DEC) ==> (xi, eta)#################################
	# use the median position of fitted star, to transform the (RA, DEC)/degree into (xi, eta)/degree
	CRA = np.median(ratmp)
 	CDEC = np.median(dectmp)
	xitmp, etatmp, status = s2t.ds2tp(ratmp, dectmp, CRA, CDEC)
	################################################################################


	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0, graErr))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0, gdecErr))*1000*3600 # #
	
 
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))# 
	rMagtmp = mr0[0]
	
	if((len(gxtmp)>0)):
	    gpmtmp = np.concatenate((gpmra, gpmraErr, gpmdec, gpmdecErr)) 
	else:
	    gpmtmp = np.array([9999, 9999, 9999, 9999])

	where_are_NaNs = np.isnan(gpmtmp)
	gpmtmp[where_are_NaNs] = 9999
	#print "gpmtmp", gpmtmp
				
	if((gpmtmp[0]<9998) or (gpmtmp[2]<9998)):
	    ndim, nwalkers = 4, 16
	    Nmcmc = 400
	    p0=np.zeros((nwalkers,ndim))
	    if(gpmtmp[0]<9998):
	    	p, cov = curve_fit(func2, xtmp*gpmtmp[0]/365.24, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    	slpxi = p[0]
	    	slpxi_sigma = np.sqrt(np.diag(cov))[0]
	    	p0[:,0] = (np.random.rand(nwalkers)*6*gpmtmp[1] + gpmtmp[0] - 3*gpmtmp[1])/365.24
	    	p0[:,1] = np.random.rand(nwalkers)*6*slpxi_sigma + slpxi - 3*slpxi_sigma
	    else:
	    	p, cov = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    	muxi = p[0]
	    	slpxi = p[1]
	    	muxi_sigma = np.sqrt(np.diag(cov))[0]
	    	slpxi_sigma = np.sqrt(np.diag(cov))[1]
	    	p0[:,0] = np.random.rand(nwalkers)*6*muxi_sigma + muxi - 3*muxi_sigma
	    	p0[:,1] = np.random.rand(nwalkers)*6*slpxi_sigma + slpxi - 3*slpxi_sigma	
	    	Nmcmc = Nmcmc + 50	
	    	
	    if(gpmtmp[2]<9998):
	    	peta, coveta = curve_fit(func2, xtmp*gpmtmp[2]/365.24, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	    	slpeta = peta[0]
	    	slpeta_sigma = np.sqrt(np.diag(coveta))[0]
	    	p0[:,2] = (np.random.rand(nwalkers)*6*gpmtmp[3] + gpmtmp[2] - 3*gpmtmp[3])/365.24
	    	p0[:,3] = np.random.rand(nwalkers)*6*slpeta_sigma + slpeta - 3*slpeta_sigma
	    else:
	    	peta, coveta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	    	mueta = peta[0]
	    	slpeta = peta[1]
	    	mueta_sigma = np.sqrt(np.diag(coveta))[0]
	    	slpeta_sigma = np.sqrt(np.diag(coveta))[1]	    	
	    	p0[:,2] = np.random.rand(nwalkers)*6*mueta_sigma + mueta - 3*mueta_sigma
	    	p0[:,3] = np.random.rand(nwalkers)*6*slpeta_sigma + slpeta - 3*slpeta_sigma
	    	Nmcmc = Nmcmc + 50	

	    #print xtmp, ytmpXI, yErrtmpXI,ytmpETA, yErrtmpETA, gpmtmp
	    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(xtmp, ytmpXI, yErrtmpXI, ytmpETA, yErrtmpETA, gpmtmp))
	    #sampler.run_mcmc(p0, Nmcmc)
	    pos, prob, state = sampler.run_mcmc(p0, 100)
	    sampler.reset()
	    sampler.run_mcmc(pos, Nmcmc, rstate0=state)
	    samples = sampler.chain[:, :, :].reshape((-1, ndim))
	    samples[:, 0] = (samples[:, 0])*365.24
	    samples[:, 2] = (samples[:, 2])*365.24
		
	    a_mcmc, b_mcmc, c_mcmc, d_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
	    #print "map:",a_mcmc[0], b_mcmc, c_mcmc[0], d_mcmc
	    muXI_mcmc = a_mcmc[0]
	    muErrXI_mcmc = (a_mcmc[1] + a_mcmc[2])/2.0
	    #b_XI_mcmc = b_mcmc[0]
	    chi2XI_mcmc = np.sum((muXI_mcmc*xtmp/365.24 + b_mcmc[0] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-2) # we should use the reduced chi2
	    muETA_mcmc = c_mcmc[0]
	    muErrETA_mcmc = (c_mcmc[1] + c_mcmc[2])/2.0
	    #b_ETA_mcmc = d_mcmc[0]
	    chi2ETA_mcmc = np.sum((muETA_mcmc*xtmp/365.24 + d_mcmc[0] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-2)
	    #gmuXI, gmuErrXI, gmuETA, gmuErrETA = gpmra, gpmraErr, gpmdec, gpmdecErr
	    if(pflag):
	        fig = corner.corner(samples, labels=[r"$\mu_{\alpha}$", r"$b_{\alpha}$", "$\mu_{\delta}$", "$b_{\delta}$"], \
		    truths=[gpmtmp[0], gpmtmp[1], gpmtmp[2], gpmtmp[3]])
	        figname = "fig_corner" + str(chunkNo) + "_%d" %idx
	        fig.savefig(root+"figs/%s.png" %figname)

	else: 
	    muXI_mcmc, muErrXI_mcmc, chi2XI_mcmc = 9999, 9999, 9999
	    muETA_mcmc, muErrETA_mcmc, chi2ETA_mcmc = 9999, 9999, 9999


	if(True): 	    
	    num0 = len(xtmp0)  + len(gxtmp) + len(xtmp1) + len(xtmp2) #X-validate PS1 + GAIA + SDSS + 2MASS points
	    s0 = 0#len(xtmp2) #
	    xpm = np.zeros([10, num0])
	    for idxX in range(0, num0): # X means cross-validation
	    	yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] + 999999.0
	    	yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] + 999999.0
	    	pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    	pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	    	xpm[0, idxX] = pX[0]*365.25
	    	xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
	    	xpm[2, idxX] = (ytmpXI[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
	    	xpm[3, idxX] = pX[1]
	    	xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - \
			    np.delete(ytmpXI, s0+idxX))**2/np.delete(yErrtmpXI, s0+idxX)**2)/(xtmp.size-3)
	    	xpm[5, idxX] = pXeta[0]*365.25
	    	xpm[6, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
	    	xpm[7, idxX] = (ytmpETA[s0+idxX] - (pXeta[0]*xtmp[s0+idxX] + pXeta[1]))
	    	xpm[8, idxX] = pXeta[1]
	    	xpm[9, idxX] = np.sum((pXeta[0]*np.delete(xtmp, s0+idxX) + pXeta[1] - \
			    np.delete(ytmpETA, s0+idxX))**2/np.delete(yErrtmpETA, s0+idxX)**2)/(xtmp.size-3)

	    	yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] - 999999.0
	    	yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] - 999999.0
	    	
	    	
	    pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
	    muXI = pXI[0]*365.25 # fit with all the data
	    muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
	    chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-2) # we should use the reduced chi2
	    b_XI = pXI[1]
	    
	    	    
	    chi2Tmp = np.concatenate((xpm[4], np.array([chi2XI])))
	    muTmp = np.concatenate((xpm[0], np.array([muXI])))
	    slpTmp = np.concatenate((xpm[3], np.array([pXI[1]])))
	    errTmp = np.concatenate((xpm[1], np.array([muErrXI])))
	    idxO = np.argmin(chi2Tmp) 
	    muXxi = muTmp[idxO] # it should be the best one
	    slpxi = slpTmp[idxO]
	    muErrXxi = errTmp[idxO]
	    chi2XI = np.min(chi2Tmp)	    
	    
	    pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
	    muETA = pETA[0]*365.25 # fit with all the data
	    muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
	    chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-2) # we should use the reduced chi2
	    b_ETA = pETA[1]	    
	    
	    chi2Tmp = np.concatenate((xpm[9], np.array([chi2ETA])))
	    muTmp = np.concatenate((xpm[5], np.array([muETA])))
	    slpTmp = np.concatenate((xpm[8], np.array([pETA[1]])))
	    errTmp = np.concatenate((xpm[6], np.array([muErrETA])))
	    idxO = np.argmin(chi2Tmp) 
	    muXeta = muTmp[idxO] # it should be the best
	    slpeta = slpTmp[idxO]
	    muErrXeta = errTmp[idxO]
	    chi2ETA = np.min(chi2Tmp)
	   

	print idx, muXI, muErrXI
	if(pflag):
	    mu = muXI
	    muErr = muErrXI
	    b =b_XI
	    yErrtmp = yErrtmpXI
	    ytmp = ytmpXI
	    muX, muErrX = gpmtmp[0], gpmtmp[1]
      
	    fig = plt.figure(figsize=(6.0,4.0))
	    ax = fig.add_subplot(111)
	    plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
	    ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
	    ax.plot(xtmp, xtmp*(mu/365.25)+b, 'r')
	    if(muX<9998): ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')
	    
	    #ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
	    #ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
	    ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
	    ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
	    ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
	    #ax.set_ylabel('$\delta-min(\delta)/mas$')
	    ax.set_xlabel('$MJD$')
	    ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
        	    %np.round(rMagtmp,decimals=2), fontsize=10)  	    
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
        	    %muErr + ",  $\mu_{g} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
	        ",  $flag = %d$" %(flag), fontsize=10)
	    ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
	    figname = "fig_PMmcmcfitting_ra" + str(chunkNo) + "_%d" %idx
	    plt.savefig(root+"figs/%s.png" %figname)
	    #plt.show()		
      
      
	    mu = muETA
	    muErr = muErrETA
	    b = b_ETA
	    yErrtmp = yErrtmpETA
	    ytmp = ytmpETA
	    muX, muErrX = gpmtmp[2], gpmtmp[3]
      
	    fig = plt.figure(figsize=(6.0,4.0))
	    ax = fig.add_subplot(111)
	    plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
	    ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
	    ax.plot(xtmp, xtmp*(mu/365.25)+b, 'r')
	    if(muX<9998): ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')
	    
	    #ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
	    #ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
	    ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
	    ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
	    #ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
	    ax.set_ylabel('$\delta-min(\delta)/mas$')
	    ax.set_xlabel('$MJD$')
	    ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
        	    %np.round(rMagtmp,decimals=2), fontsize=10)
	    ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
        	    %muErr + ",  $\mu_{g} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
	        ",  $flag = %d$" %(flag), fontsize=10)
	    ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
	    figname = "fig_PMmcmcfitting_dec" + str(chunkNo) + "_%d" %idx
	    plt.savefig(root+"figs/%s.png" %figname)
	    #plt.show()				
		
        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        #ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
	#print idx
        return uniqueID, ra0[0], dec0[0], rMagtmp, muXI, muErrXI, muXxi, muErrXxi, chi2XI, \
	    muXI_mcmc, muErrXI_mcmc, chi2XI_mcmc, gpmtmp[0], gpmtmp[1], \
	    muETA, muErrETA, muXeta, muErrXeta, chi2ETA, muETA_mcmc, muErrETA_mcmc, \
	    chi2ETA_mcmc, gpmtmp[2], gpmtmp[3], len(xtmp0), len(xtmp), flag
    else:
	return idx, uniqueID

##############fitting PM for each star parallelly, + GAIA DR2, with MCMC ################
def fittingPM_mcmc2(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
        obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
        dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
        gpmra, gpmdec, gpmraErr, gpmdecErr, \
        pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
    ind0 = obj_id0 == uniqueID
    xtmp0 = mjd0[ind0] 
    if(len(xtmp0)>2):
	ra0 = ra0[ind0]
	raErr0 = raErr0[ind0]
	dec0 = dec0[ind0]
	decErr0 = decErr0[ind0]
	mr0 = mr0[ind0]
	#print uniqueID, ra0, dec0, mr0


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
	gpmraErr = gpmraErr[gind]
	gpmdecErr = gpmdecErr[gind]
	if(len(gxtmp)>0): flag = flag + 20

	ratmp = np.concatenate((ra1, ra2, ra0, gra))
	dectmp = np.concatenate((dec1, dec2, dec0, gdec))
	#print ratmp, flag

	#########################(RA, DEC) ==> (xi, eta)#################################
	# use the median position of fitted star, to transform the (RA, DEC)/degree into (xi, eta)/degree
	CRA = np.median(ratmp)
 	CDEC = np.median(dectmp)
	xitmp, etatmp, status = s2t.ds2tp(ratmp, dectmp, CRA, CDEC)
	################################################################################


	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0, graErr))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0, gdecErr))*1000*3600 # #
	
 
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))# 
	rMagtmp = mr0[0]
	
	if((len(gxtmp)>0)):
	    gpmtmp = np.concatenate((gpmra, gpmraErr, gpmdec, gpmdecErr)) 
	else:
	    gpmtmp = np.array([9999, 0, 9999, 0])

	where_are_NaNs = np.isnan(gpmtmp)
	gpmtmp[where_are_NaNs] = 9999
	#print "gpmtmp", gpmtmp
				
	if((len(xtmp)>2)):
		ndim, nwalkers = 4, 30
		Nmcmc = 400
		p0=np.zeros((nwalkers,ndim))
		if((gpmtmp[0]<9998) or (gpmtmp[2]<9998)):
			p0[:,0] = (np.random.rand(nwalkers)*2*gpmtmp[1] + gpmtmp[0] - gpmtmp[1])/365.24
			p0[:,2] = (np.random.rand(nwalkers)*2*gpmtmp[3] + gpmtmp[2] - gpmtmp[3])/365.24
		else:
			p0[:,0] = np.random.rand(nwalkers)*4-2
			p0[:,2] = np.random.rand(nwalkers)*4-2
			Nmcmc = 1000
		
		p0[:,1] = np.random.rand(nwalkers)*100-50
		p0[:,3] = np.random.rand(nwalkers)*100-50
		

		#print xtmp, ytmpXI, yErrtmpXI,ytmpETA, yErrtmpETA, gpmtmp
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(xtmp, ytmpXI, yErrtmpXI, ytmpETA, yErrtmpETA, gpmtmp))
		sampler.run_mcmc(p0, Nmcmc)
		samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
		samples[:, 0] = (samples[:, 0])*365.24
		samples[:, 2] = (samples[:, 2])*365.24
		
			
		a_mcmc, b_mcmc, c_mcmc, d_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
		#print "map:",a_mcmc[0], b_mcmc, c_mcmc[0], d_mcmc
		muXI = a_mcmc[0]
		muErrXI = (a_mcmc[1] + a_mcmc[2])/2.0
		chi2XI = np.sum((muXI*xtmp/365.24 + b_mcmc[0] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-2) # we should use the reduced chi2
		muETA = c_mcmc[0]
		muErrETA = (c_mcmc[1] + c_mcmc[2])/2.0
		chi2ETA = np.sum((muETA*xtmp/365.24 + d_mcmc[0] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-2)
		#gmuXI, gmuErrXI, gmuETA, gmuErrETA = gpmra, gpmraErr, gpmdec, gpmdecErr
		print muXI, chi2XI, muETA, chi2ETA
		if(idx>10): pflag=False
		if(pflag):
			fig = corner.corner(samples, labels=["$\mu_{\alpha}$", "$b_{\alpha}$", "$\mu_{\delta}$", "$b_{\delta}$"], \
			    truths=[gpmtmp[0], gpmtmp[1], gpmtmp[2], gpmtmp[3]])
			figname = "fig_corner" + str(chunkNo) + "_%d" %idx
			fig.savefig(root+"figs/%s.png" %figname)

			mu = a_mcmc[0]
			muErr = (a_mcmc[1] + a_mcmc[2])/2.0
			yErrtmp = yErrtmpXI
			ytmp = ytmpXI
			muX, muErrX = gpmtmp[0], gpmtmp[1]

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(mu/365.25)+b_mcmc[0], 'r')
			if(muX<9998): ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')
			
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			#ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)  	    
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{g} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMmcmcfitting_ra" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()		

		
			mu = c_mcmc[0]
			muErr = (c_mcmc[1] + c_mcmc[2])/2.0
			yErrtmp = yErrtmpETA
			ytmp = ytmpETA
			muX, muErrX = gpmtmp[2], gpmtmp[3]

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(mu/365.25)+d_mcmc[0], 'r')
			if(muX<9998): ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')
			
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{g} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMmcmcfitting_dec" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/%s.png" %figname)
			#plt.show()			
		
        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        #ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
	#print idx
        return uniqueID, ra0[0], dec0[0], rMagtmp, muXI, muErrXI, chi2XI, gpmra, gpmraErr, \
	    muETA, muErrETA, chi2ETA, gpmdec, gpmdecErr, len(xtmp0), len(xtmp), flag
    else:
	return idx, uniqueID
##############fitting PM for each star parallelly, + GAIA, speed up basedon fittingPM4 ################
def fittingPM4QSOs(packParameterList):
    idx, uniqueID, obj_id0, ra0, raErr0, dec0, decErr0, mjd0, mr0, \
	obj_id1, ra1, raErr1, dec1, decErr1, mjd1, obj_id2, ra2, raErr2, \
	dec2, decErr2, mjd2, gobj_id, gra, graErr, gdec, gdecErr, gmjd, \
	pobj_id, pmrapv3, pmdecpv3, epmrapv3, epmdecpv3, pflag, chunkNo, CRA, CDEC, root = packParameterList
    #t0 = time()
    flag = 0  #record the PS1 
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

	pind = pobj_id == uniqueID



	xitmp = np.concatenate((ra1, ra2, ra0, gra))#   #
	etatmp = np.concatenate((dec1, dec2, dec0, gdec))#  # 

	ytmpXI = (xitmp-np.min(xitmp))*1000*3600
	yErrtmpXI = np.concatenate((raErr1, \
	    raErr2, raErr0, graErr))*1000*3600 # # 

	ytmpETA = (etatmp-np.min(etatmp))*1000*3600
	yErrtmpETA = np.concatenate((decErr1, \
	    decErr2, decErr0, gdecErr))*1000*3600 # #
	
	xtmp = np.concatenate((xtmp1,xtmp2, xtmp0, gxtmp))#  

	#print "preparing time:", time()-t0
	if((len(xtmp)>2)):
		rMagtmp = mr0[0]
		num0 = len(xtmp0)  + len(gxtmp) #X-validate PS1 + GAIA points
		s0 = len(xtmp1) + len(xtmp2) #	
		#print "starIdx, num0, chunkNo:", idx, num0,chunkNo
		xpm = np.zeros([10, num0])
		for idxX in range(0, num0): # X means cross-validation
			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] + 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] + 999999.0
			pX, covX = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
			pXeta, covXeta = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
			xpm[0, idxX] = pX[0]*365.25
			xpm[1, idxX] = np.sqrt(np.diag(covX))[0]*365.25
			xpm[2, idxX] = (ytmpXI[s0+idxX] - (pX[0]*xtmp[s0+idxX] + pX[1]))
			xpm[3, idxX] = pX[1]
			xpm[4, idxX] = np.sum((pX[0]*np.delete(xtmp, s0+idxX) + pX[1] - \
			    np.delete(ytmpXI, s0+idxX))**2/np.delete(yErrtmpXI, s0+idxX)**2)/(xtmp.size-2)

			xpm[5, idxX] = pXeta[0]*365.25
			xpm[6, idxX] = np.sqrt(np.diag(covXeta))[0]*365.25
			xpm[7, idxX] = (ytmpETA[s0+idxX] - (pXeta[0]*xtmp[s0+idxX] + pXeta[1]))
			xpm[8, idxX] = pXeta[1]
			xpm[9, idxX] = np.sum((pXeta[0]*np.delete(xtmp, s0+idxX) + pXeta[1] - \
			    np.delete(ytmpETA, s0+idxX))**2/np.delete(yErrtmpETA, s0+idxX)**2)/(xtmp.size-2)

			yErrtmpXI[s0+idxX] = yErrtmpXI[s0+idxX] - 999999.0
			yErrtmpETA[s0+idxX] = yErrtmpETA[s0+idxX] - 999999.0

		pXI, covXI = curve_fit(func, xtmp, ytmpXI, sigma=yErrtmpXI, absolute_sigma=True)
		muXI = pXI[0]*365.25 # fit with all the data
		muErrXI = np.sqrt(np.diag(covXI))[0]*365.25
		chi2XI = np.sum((pXI[0]*xtmp + pXI[1] - ytmpXI)**2/yErrtmpXI**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[4], np.array([chi2XI])))
	        muTmp = np.concatenate((xpm[0], np.array([muXI])))
	        slpTmp = np.concatenate((xpm[3], np.array([pXI[1]])))
	        errTmp = np.concatenate((xpm[1], np.array([muErrXI])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXxi = muTmp[idxO] # it should be the best one
	        slpxi = slpTmp[idxO]
		muErrXxi = errTmp[idxO]
	        chi2XI = np.min(chi2Tmp)
		
		muXIog = xpm[0][num0-1] # the last one should be the fit without GAIA
		muErrXIog = xpm[1][num0-1]

		pETA, covETA = curve_fit(func, xtmp, ytmpETA, sigma=yErrtmpETA, absolute_sigma=True)
		muETA = pETA[0]*365.25 # fit with all the data
		muErrETA = np.sqrt(np.diag(covETA))[0]*365.25
		chi2ETA = np.sum((pETA[0]*xtmp + pETA[1] - ytmpETA)**2/yErrtmpETA**2)/(xtmp.size-1) # we should use the reduced chi2

	        chi2Tmp = np.concatenate((xpm[9], np.array([chi2ETA])))
	        muTmp = np.concatenate((xpm[5], np.array([muETA])))
	        slpTmp = np.concatenate((xpm[8], np.array([pETA[1]])))
	        errTmp = np.concatenate((xpm[6], np.array([muErrETA])))
 	        idxO = np.argmin(chi2Tmp) 
	        muXeta = muTmp[idxO] # it should be the best
	        slpeta = slpTmp[idxO]
		muErrXeta = errTmp[idxO]
	        chi2ETA = np.min(chi2Tmp)

		muETAog = xpm[5][num0-1] # the last one should be the fit without GAIA
		muErrETAog = xpm[6][num0-1]

		if(pflag):
			p = pETA#pXI
			mu = muETA#muXI
			muErr = muErrETA#muErrXI
			yErrtmp = yErrtmpETA#yErrtmpXI
			ytmp = ytmpETA#ytmpXI
			muX = muXeta#muXxi
			muErrX = muErrXeta#muErrXxi
			#R2 = R2ETA#R2XI
			chi2 = chi2ETA#chi2XI

			fig = plt.figure(figsize=(6.0,4.0))
			ax = fig.add_subplot(111)
    			plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
			ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
			ax.plot(xtmp, xtmp*(muX/365.25)+(np.median(ytmp)-muX/365.25*np.median(xtmp)), 'b--')
			ax.plot(xtmp, xtmp*(pmdecpv3[pind]/365.25)+(np.median(ytmp)-pmdecpv3[pind]/365.25*np.median(xtmp)), 'k--')

			ax.text(np.min(xtmp)+50,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=10)
			#ax.text(np.min(xtmp)+50,np.max(ytmp)*0.82,'$R^2 = %0.2f$'% R2, fontsize=10)
			ax.plot(xtmp, xtmp*p[0]+p[1], 'r-')
			ax.set_ylim([int(-1*np.max(ytmp)/4), np.max(ytmp)*1.4])	
			ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
			#ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr)$')
			ax.set_ylabel('$\delta-min(\delta)/mas$')
			ax.set_xlabel('$MJD$')
			ax.set_xlim(np.min(xtmp)-100,np.max(xtmp)+100) 
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.25, "$objID = %d$" %uniqueID, fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.15, "$rmag = %0.2f$" \
		    	    %np.round(rMagtmp,decimals=2), fontsize=10)
			ax.text(np.min(xtmp)+50, np.max(ytmp)*1.05, "$\mu = %0.2f$" %(mu)+"$\pm%0.2f$" \
		    	    %muErr + ",  $\mu_{X} = %0.2f$" %(muX)+"$\pm%0.2f$"%muErrX + \
			    ",  $\mu_{pv3} = %0.2f$" %(pmdecpv3[pind])+"$\pm%0.2f$"%epmdecpv3[pind] + \
			    ",  $flag = %d$" %(flag), fontsize=10)
			ax.set_title('Proper Motions with Season-AVG PS1 and Other Surveys', fontsize=12)
			figname = "fig_PMfitting_ETA" + str(chunkNo) + "_%d" %idx
			plt.savefig(root+"figs/QSOs/%s.png" %figname)
			#plt.show()

        #########################(xi, eta) ==> (RA, DEC)#################################
        # transform the (xi, eta)/radians into (RA, DEC)/degree
        ra, dec = s2t.dtp2s(np.radians(ra0[0]), np.radians(dec0[0]), CRA, CDEC)
        #################################################################################
        return uniqueID, ra, dec, rMagtmp, muXI, muErrXI, muXxi, muErrXxi, muXIog, muErrXIog, chi2XI, \
	    muETA, muErrETA, muXeta, muErrXeta, muETAog, muErrETAog, chi2ETA, len(xtmp), flag
    else:
	return idx, uniqueID

##############fitting PM in different region################
def fitPMdist(obj_id, ra,dec, ra0, dec0, decErr, mjd, mr, root):
    raTmp, decTmp = s2t.dtp2s(np.radians(ra), np.radians(dec), ra0, dec0)
    ra1 = np.min(raTmp) + (np.max(raTmp)-np.min(raTmp))/2 
    dec1 = np.min(decTmp) + (np.max(decTmp)-np.min(decTmp))/2
    distTmp = sphdist(ra1, dec1, raTmp, decTmp)
    distMask = np.argsort(distTmp)
    #distMask = distMask[0:100]
    obj_id = obj_id[distMask]
    ra = ra[distMask]
    dec = dec[distMask]
    #raErr = raErr[distMask]
    decErr = decErr[distMask]
    mr = mr[distMask]
    mjd = mjd[distMask]
    uniqueStar = np.unique(obj_id)
    print "len(uniqueStar)", len(uniqueStar)
    print "uniqueStar[0:20]", uniqueStar[0:30]
    #uniqueStar2 = np.array([-6777811886076157207, -6777811886076157083, -6777811886076157069, -6777811886076048565, \
    #	-6759797487566675242, -6759797487566675212, -6759797487566675202, -6759797487566675165, -6759797487566675148, -6759797487566675144])
    '''uniqueStar2 = np.array([-6922050219454086081, \
 	-6904053413130931635, -6904053413130931533, -6904053413130931524, \
 	-6904053413130931464, -6904053413130931450, \
 	-6904053413130931408, -6904053413130931376, -6904053413130931358, \
 	-6904053413130931306, -6886039014621449656, -6886039014621449641, \
 	-6886039014621449582, -6886039014621449575, -6886039014621449544, \
 	-6886039014621449537, -6886039014621449532, -6886039014621449516, \
 	-6886039014621449503, -6886039014621449473, -6886039014621449470, \
 	-6886039014621449457, -6886039014621449386])
    uniqueStar2 = np.array([-6814086973699879558, -6814086973699879541, -6814086973699879475, -6814086973699879453, \
	-6814086973699879367, -6814086973699879315, -6814086973699879273, -6814086973699879203, \
	-6814086973699879160, -6814086973699657457, -6814086973699657416, -6814086973699613232])
    uniqueStar2 = np.array([-4799834448172286076,-4799834448172286072,-4799834448172285748 \
	,-4799834448172285514,-4799834448172285441,-4799834448172285307 \
	,-4799834448172284949,-4799834448172284823 \
	,-4799834448172284526,-4799834448172249305,-4799834448172249055 \
	,-4799834448172248849,-4799834448171775835,-4799834448171774940 \
	,-4799834448171774391,-4799834448171774084,-4799834448171742013 \
	,-4799834448171633874,-4799834448171611781,-4799834448171427317 \
	,-4799834448171157104,-4799834448171033515,-4799834448171031649 \
	,-4799834448171018557,-4799834448170954263,-4799834448170954237 \
	,-4799834448170954223,-4799834448170014230])'''
    uniqueStar2 = np.array([-4871962410954520621,-4871962410954520617 \
	,-4871962410954520549,-4871962410954520532,-4871962410954520531 \
	,-4871962410954520433,-4871962410954520372,-4871962410954520368 \
	,-4871962410954519821,-4871962410954519484,-4871962410954519338 \
	,-4871962410954519037,-4871962410954470165,-4871962410954470116 \
	,-4871962410954469682,-4871962410954469637,-4871962410954469525 \
	,-4871962410954228152,-4871962410954209854,-4871962410954190579 \
	,-4871962410954190482,-4871962410954172430,-4871962410954172429 \
	,-4871962410953931840,-4871962410953931703,-4871962410953861782 \
	,-4871962410953861340,-4871962410953134470])

    R2 = []
    idx=0
    count = 0
    fig, axs = plt.subplots(nrows=5, ncols=2, sharex=True, figsize=(8.0,12.0))
    #plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
    plt.subplots_adjust(left=0.07, bottom=0.05, right=0.99, top=0.99, wspace=0.16, hspace=0.)
    while (count<10):
	ind = obj_id == uniqueStar2[idx]
	print np.sum(ind)
	xtmp = mjd[ind]
	tmp = ra[ind]
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = decErr[ind]*1000*3600
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	if((len(xtmp)>2)):
		ax = axs.ravel()[count]
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')
		ax.set_xlim(55400, 57000)#(55000,56500) 
		#x_extra = np.concatenate((xtmp,xtmp[-1:]))
		#y_extra = np.concatenate((ytmp, ytmp[-1:]))
		#weights = [1 for i in range(xtmp.size)]
		#weights.append(sys.float_info.epsilon)
		#p,cov = np.polyfit(xtmp, ytmp, 1, w=yErrtmp, cov=True)
		p, cov = curve_fit(func, xtmp, ytmp, sigma=yErrtmp, absolute_sigma=True)
		# coefficient of determination, plot text
		variance = np.var(ytmp)
		#residuals = np.var([(p[0]*xx + p[1] - yy)  for xx,yy in zip(xtmp,ytmp)])
		residuals = np.var(ytmp - p[0]*xtmp + p[1])
		Rtmp = residuals/variance		
		Rsqr = np.round(np.abs(1-Rtmp), decimals=2) 
		R2.append(Rsqr)
		chi2 = np.sum((p[0]*xtmp + p[1] - ytmp)**2/yErrtmp**2)#/(xtmp.size-1)
		ax.text(55600,np.max(ytmp)*0.95,'$\chi^2 = %0.2f$'% chi2, fontsize=8)
		ax.text(55600,np.max(ytmp)*0.85,'$R^2 = %0.2f$'% Rsqr, fontsize=8)
		ax.plot(xtmp, xtmp*p[0]+p[1], '-')
		ax.set_ylim([int(-1*np.max(ytmp)/3), np.max(ytmp)*1.4])	
		ax.set_yticks([0, int(np.max(ytmp)/2), int(np.max(ytmp))])
		ax.set_ylabel('$\eta-min(\eta)/mas$')
		if(count>7):	ax.set_xlabel('$MJD')
		
		ax.text(55600, np.max(ytmp)*1.25, "$objID=%d$" %objIDtmp, fontsize=8)
		ax.text(55600, np.max(ytmp)*1.15, "$rmag=%0.2f$" %np.round(rMagtmp,decimals=2), fontsize=8)
		ax.text(55600, np.max(ytmp)*1.05, "$\mu_{\eta}=%0.2f$" %(p[0]*360)+"$\pm$%0.2f" \
		    %(np.sqrt(np.diag(cov))[0]*360), fontsize=8)
		#ax.set_title('Postional Variation of Star in Different Epochs', fontsize=12)
		count = count + 1
	idx = idx + 1
    figname = "fig_" + "%d" %ra1+ "_%d" %dec1
    plt.savefig(root+"figs/%s_pv3_uncorr.png" %figname)

def scatterMagMu(magPS1, magSDSS, dltMuPS1, dltMuSDSS, minr, maxr, dltr, lg, txt, root):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(5.5,5.0))
    plt.subplots_adjust(left=0.135, bottom=0.12, right=0.99, top=0.99, wspace=0.16, hspace=0.)
    #####log######
    if(lg):
        dltMuPS1 = np.log10(dltMuPS1)
        dltMuSDSS = np.log10(dltMuSDSS)
    magBreakAt = np.arange(minr, maxr, dltr)
    #magBreakAt = np.arange(15.2, 15.8, 0.2)
    medianMag = np.zeros(len(magBreakAt)+1)
    medianDeltaMu = np.zeros(len(magBreakAt)+1)
    medianMagSDSS = np.zeros(len(magBreakAt)+1)
    medianDeltaMuSDSS = np.zeros(len(magBreakAt)+1)
    errDeltaMuPS1 = np.zeros(len(magBreakAt)+1)
    errDeltaMuSDSS = np.zeros(len(magBreakAt)+1)
    for idx in range(0, len(magBreakAt)+1):
	if(idx==0):
	    var = magBreakAt[idx]
	    magInBin = magPS1 < var
	    magInBinSDSS = magSDSS < var
	elif(idx==len(magBreakAt)):
	    var = magBreakAt[idx-1]
	    magInBin = magPS1 >= var
	    magInBinSDSS = magSDSS >= var
	else:
	    var = magBreakAt[idx]
	    magInBin = (magPS1 >= magBreakAt[idx-1]) & (magPS1 < var)
            magInBinSDSS = (magSDSS >= magBreakAt[idx-1]) & (magSDSS < var)

	print "idx, sum(magInBin):", idx, sum(magInBin), sum(magInBinSDSS)
	medianMag[idx] = np.median(magPS1[magInBin])
	medianMagSDSS[idx] = np.median(magSDSS[magInBinSDSS])
	medianDeltaMu[idx] = np.median(dltMuPS1[magInBin])
	medianDeltaMuSDSS[idx] = np.median(dltMuSDSS[magInBinSDSS])
 	errDeltaMuPS1[idx] = (0.741*(np.percentile(dltMuPS1[magInBin], 75) - np.percentile(dltMuPS1[magInBin], 25)))
 	errDeltaMuSDSS[idx] = (0.741*(np.percentile(dltMuSDSS[magInBinSDSS], 75) - np.percentile(dltMuSDSS[magInBinSDSS], 25)))

        print np.sum(magInBin), np.sum(magInBinSDSS)


    ax = axs.ravel()[0]
    mTmp = np.median(dltMuPS1[(magPS1<18)&(magPS1>14)])
    print "the mean value of muRaERR (mas/yr):", mTmp
    if(lg):
	#if(txt=="GPS1"): ymax = 0.95
	#else: ymax = 0.95
	ymax = 1.24
        ax.plot([0, 25], [mTmp, mTmp], 'k--', linewidth=1.5)
        ax.text(15, ymax*0.8, "$" + txt + "\ Stars$")#"ONLY PS1 "+txt)# 
        ax.set_ylim([0,ymax])
        #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
        ax.set_yticks(np.arange(0,ymax,0.25))
        ax.set_ylabel(r'$\log_{10}(\Delta \mu_{\alpha*}(mas/yr))$', fontsize=18)
    else:
	ymax = 11.0
	ax.set_yscale('log')
        ax.plot([0, 25], [mTmp, mTmp], 'k--', linewidth=1.5)
        ax.text(15.0, ymax*0.6, "$" + txt + "\ Stars$")# "OUR PM(PS1 no zero-point shifts)")#+txt
        ax.set_ylim([0.5,ymax])
        ax.set_yticks([1.0, 2.0, 4.0, 8.0])#ax.set_yticks(np.arange(0,ymax,2.0)) #0.5
	ax.get_yaxis().set_major_formatter(ScalarFormatter())
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.01f'))#'%.02f'  '%d'
	#ax.set_yticklabels(["2", "4", "6", "8", "10"])
        #ax.set_ylabel(r'$\Delta \mu_{\alpha}(mas/yr)$')
	#ax.set_ylabel(r'$\boldsymbol{\delta}\mu_{\alpha*}(mas/yr)$', fontsize='large')	
	ax.set_ylabel(r'$\epsilon_{\mu_{\alpha*}}(mas/yr)$', fontsize=18)

    ax.scatter(magPS1, dltMuPS1, color='b', marker='.', lw=2.0, alpha=0.02, rasterized=True) #0.005
    ax.errorbar(medianMag, medianDeltaMu, errDeltaMuPS1, color = 'r', elinewidth=2)    
    ax.minorticks_on()
    ax.text(13.5, 0.7, '$Dashed:\ ' + str('%.2f'%(mTmp)) + '\ mas/yr$', fontsize=12)
    ax.set_xlim([13.2, 20.1])#([14, 19.5])#
    ax.set_xlabel('')


    ax = axs.ravel()[1]
    mTmp = np.median(dltMuSDSS[(magSDSS<18)&(magSDSS>14)])
    print "the mean value of muDecERR (mas/yr):", mTmp
    if(lg):
	#if(txt=="GPS1"): ymax = 0.95
	#else: ymax = 0.95
	ymax = 1.24
        ax.plot([0, 25], [mTmp, mTmp], 'k--', linewidth=1.5)
        ax.text(15, ymax*0.8, "$" + txt + "\ Stars$")#"PS1 + SDSS "+txt)  #
        ax.set_ylim([0,ymax])
        #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
        ax.set_yticks(np.arange(0,ymax,0.25))#5
        ax.set_ylabel(r'$\log_{10}(\Delta \mu_{\delta}(mas/yr))$', fontsize=18)
    else:
	ymax = 11
	ax.set_yscale('log')
        ax.plot([0, 25], [mTmp, mTmp], 'k--', linewidth=1.5)
        ax.text(15.0, ymax*0.6, "$" + txt + "\ Stars$")#"PS1 + SDSS "+txt) 
        ax.set_ylim([0.5,ymax])
        ax.set_yticks([1.0, 2.0, 4.0, 8])#(np.arange(1,ymax,2.0)) #0.5
	ax.get_yaxis().set_major_formatter(ScalarFormatter())
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.01f')) #'%.02f'
        #ax.set_ylabel('$\Delta \mu_{\delta}(mas/yr)$')
	ax.set_ylabel(r'$\epsilon_{\mu_{\delta}}(mas/yr)$', fontsize=18)

    ax.scatter(magSDSS, dltMuSDSS, color='b', marker='.', lw=2.0, alpha=0.02, rasterized=True) #0.01
    ax.errorbar(medianMagSDSS, medianDeltaMuSDSS, errDeltaMuSDSS, color='r', elinewidth=2)
    ax.minorticks_on()	
    ax.text(13.5, 0.7, '$Dashed:\ ' + str('%.2f'%(mTmp)) + '\ mas/yr$', fontsize=12)
    ax.set_xlim([13.2, 20.1])##([minr-dltr+0.2,maxr+dltr])
    ax.set_xlabel('$r(mag)$', fontsize=18)
    plt.savefig(root+"figs/rmag_muErr_"+txt+"_mrcal2.pdf")


def scatterMagMu2(mag1, mag2, mag3, dltMu1, dltMu2,dltMu3, lg, txt, root):
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(5.5,8.0))
    plt.subplots_adjust(left=0.135, bottom=0.08, right=0.99, top=0.99, wspace=0.16, hspace=0.)
    #####log######
    if(lg):
        dltMu1 = np.log10(dltMu1)
        dltMu2 = np.log10(dltMu2)
        dltMu3 = np.log10(dltMu3)
    magBreakAt = np.arange(15.2, 16.8, 0.2)
    #magBreakAt = np.arange(15.2, 15.8, 0.2)
    medianMag = np.zeros(len(magBreakAt)+1)
    medianDeltaMu = np.zeros(len(magBreakAt)+1)
    medianMag2 = np.zeros(len(magBreakAt)+1)
    medianDeltaMu2 = np.zeros(len(magBreakAt)+1)
    medianMag3 = np.zeros(len(magBreakAt)+1)
    medianDeltaMu3 = np.zeros(len(magBreakAt)+1)
    errDeltaMu1 = np.zeros(len(magBreakAt)+1)
    errDeltaMu2 = np.zeros(len(magBreakAt)+1)
    errDeltaMu3 = np.zeros(len(magBreakAt)+1)
    for idx in range(0, len(magBreakAt)+1):
	if(idx==0):
	    var = magBreakAt[idx]
	    magInBin = mag1 < var
	    magInBin2 = mag2 < var
	    magInBin3 = mag3 < var
	elif(idx==len(magBreakAt)):
	    var = magBreakAt[idx-1]
	    magInBin = mag1 >= var
	    magInBin2 = mag2 >= var
	    magInBin3 = mag3 >= var
	else:
	    var = magBreakAt[idx]
	    magInBin = (mag1 >= magBreakAt[idx-1]) & (mag1 < var)
            magInBin2 = (mag2 >= magBreakAt[idx-1]) & (mag2 < var)
            magInBin3 = (mag3 >= magBreakAt[idx-1]) & (mag3 < var)

	medianMag[idx] = np.median(mag1[magInBin])
	medianMag2[idx] = np.median(mag2[magInBin2])
	medianMag3[idx] = np.median(mag3[magInBin3])
	medianDeltaMu[idx] = np.median(dltMu1[magInBin])
	medianDeltaMu2[idx] = np.median(dltMu2[magInBin2])
	medianDeltaMu3[idx] = np.median(dltMu3[magInBin3])
 	errDeltaMu1[idx] = (0.741*(np.percentile(dltMu1[magInBin], 75) - np.percentile(dltMu1[magInBin], 25)))
 	errDeltaMu2[idx] = (0.741*(np.percentile(dltMu2[magInBin2], 75) - np.percentile(dltMu2[magInBin2], 25)))
 	errDeltaMu3[idx] = (0.741*(np.percentile(dltMu3[magInBin3], 75) - np.percentile(dltMu3[magInBin3], 25)))

        print np.sum(magInBin), np.sum(magInBin2)
    ax = axs.ravel()[0]
    ax.scatter(mag1, dltMu1, color='b', marker='.', lw=2.0, alpha=0.1) #0.01
    ax.errorbar(medianMag, medianDeltaMu, errDeltaMu1, color = 'r', elinewidth=1)
    if(lg):
        ax.plot([0, 25], [0.18, 0.18], 'k--', linewidth=1.5)
        ax.text(15.0, 1.25, "OUR PM(PS1+SDSS)")#"ONLY PS1 "+txt)# 
        ax.set_ylim([0,1.5])
        #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
        ax.set_yticks(np.arange(0,1.4,0.25))
        ax.set_ylabel('$\log_{10}(\Delta \mu_{\delta}(mas/yr))$')
    else:
        ax.plot([0, 25], [1.5, 1.5], 'k--', linewidth=1.5)
        ax.text(15.3, 2.5, "OUR PM(PS1)")# "OURs PM(PS1 no zero-point shifts)")#+txt
        ax.set_ylim([0.5,3])
        ax.set_yticks(np.arange(0.5,2.9,0.5))
        ax.set_ylabel('$\Delta \mu_{\delta}(mas/yr)$')	
    ax.set_xlim([13,23])#([14, 19.5])#
    ax.set_xlabel('')

    ax = axs.ravel()[1]
    ax.scatter(mag2, dltMu2, color='b', marker='.', lw=2.0, alpha=0.1) #0.01
    ax.errorbar(medianMag2, medianDeltaMu2, errDeltaMu2, color='r', elinewidth=1)

    if(lg):
        ax.plot([0, 25], [0.18, 0.18], 'k--', linewidth=1.5)
        ax.text(15.0, 1.25, "PS1's PM (Gene)")#"PS1 + SDSS "+txt)  #
        ax.set_ylim([0,1.5])
        #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
        ax.set_yticks(np.arange(0,1.4,0.25))
        ax.set_ylabel('$\log_{10}(\Delta \mu_{\delta}(mas/yr))$')
    else:
        ax.plot([0, 25], [1.5, 1.5], 'k--', linewidth=1.5)
        ax.text(15.3, 2.5, "PS1's PM (Gene)")#"PS1 + SDSS "+txt) 
        ax.set_ylim([0.5,3])
        ax.set_yticks(np.arange(0.5,2.9,0.5))
        ax.set_ylabel('$\Delta \mu_{\delta}(mas/yr)$')
    ax.set_xlim([13,23])#([14, 19.5])#
    ax.set_xlabel('')	

    ax = axs.ravel()[2]
    ax.scatter(mag3, dltMu3, color='b', marker='.', lw=2.0, alpha=0.1) #0.01
    ax.errorbar(medianMag3, medianDeltaMu3, errDeltaMu3, color='r', elinewidth=1)

    if(lg):
        ax.plot([0, 25], [0.18, 0.18], 'k--', linewidth=1.5)
        ax.text(15.0, 1.25, "Pal5 PM (Tobias)")#"PS1 + SDSS "+txt)  #
        ax.set_ylim([0,1.5])
        #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
        ax.set_yticks(np.arange(0,1.4,0.25))
        ax.set_ylabel('$\log_{10}(\Delta \mu_{\delta}(mas/yr))$')
    else:
        ax.plot([0, 25], [1.5, 1.5], 'k--', linewidth=1.5)
        ax.text(15.3, 2.5, "Pal5 PM (Tobias)")#"PS1 + SDSS "+txt) 
        ax.set_ylim([0.5,3])
        ax.set_yticks(np.arange(0.5,2.9,0.5))
        ax.set_ylabel('$\Delta \mu_{\delta}(mas/yr)$')

    ax.set_xlim([14.8,17.1])#([14, 19.5])#
    ax.set_xlabel('$r(mag)$')
    plt.savefig(root+"figs/rmag_muErr_OursPS1_Gene_Tobias_"+txt+".png")


def scatterMagMu3(mr, mura, mudec, minr, maxr, dltr, txt, filename):

    binnedRA = useful.med_bin(mr, mura, minr, maxr, dltr)
    binnedDEC = useful.med_bin(mr, mudec, minr, maxr, dltr)
    gdraREF, gddecREF = 0, 0


    good = binnedRA[4] > 20
    print "mpmraErr, mpmdecErr:", np.mean(binnedRA[2][good]), np.mean(binnedDEC[2][good])

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(5.5,5.0))
    plt.subplots_adjust(left=0.155, bottom=0.12, right=0.99, top=0.99, wspace=0.16, hspace=0.)
    ax = axs.ravel()[0]
    ax.scatter(mr, mura, color='b', marker='.', lw=2.0, alpha=0.05, rasterized=True)
    ax.errorbar(binnedRA[0][good], binnedRA[1][good], yerr=binnedRA[2][good], fmt='-o', color='r', ecolor='r', elinewidth=2.0)
    ax.plot([0, 24], [gdraREF, gdraREF], 'k--')
    ax.minorticks_on()
    ax.set_xlim([minr, maxr+0.1])
    #ax.set_ylim([gdraREF-50*max(binnedRA[3][good]), gdraREF+50*max(binnedRA[3][good])])
    ax.set_ylim([-11., 11.])
    ax.set_xlabel('')
    ax.set_ylabel(txt[0], fontsize=18)

    ax = axs.ravel()[1]
    ax.scatter(mr, mudec, color='b', marker='.', lw=2.0, alpha=0.05, rasterized=True)
    ax.errorbar(binnedDEC[0][good], binnedDEC[1][good], yerr=binnedDEC[2][good], fmt='-o', color='r', ecolor='r', elinewidth=2.0)
    ax.plot([0, 24], [gddecREF, gddecREF], 'k--')

    ax.minorticks_on()	
    ax.set_xlim([minr, maxr+0.1])
    #ax.set_ylim([gddecREF-50*max(binnedDEC[3][good]), gddecREF+50*max(binnedDEC[3][good])])
    ax.set_ylim([-11.0, 11.0])
    #ax.text(minr+0.5, gddecREF-30*max(binnedDEC[3][good]), txt[2])
    ax.text(minr+0.5,-7.5, txt[2], fontsize=16)
    ax.set_xlabel('$r(mag)$', fontsize=18)
    ax.set_ylabel(txt[1], fontsize=18)
    plt.savefig(filename)


def histMu(muPS1, muSDSS, txt, mr0, mr1, fname):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(4.5,6.0))
    plt.subplots_adjust(left=0.15, bottom=0.12, right=0.98, top=0.99, wspace=0.16, hspace=0.)
    Nbin = 15
    xmin, xmax = -13, 13
    ax = axs.ravel()[0]  
    histo = ax.hist(muPS1, Nbin, [xmin, xmax],color='k')
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
    mu = np.median(muPS1)
    sigma = (0.741*(np.percentile(muPS1, 75) - np.percentile(muPS1, 25)))
    print mu, sigma
    ax.plot([0, 0], [0, 19000], '--', markersize=2.5, color='y')
    ax.set_ylim([0, np.max(histo[0])*1.1])	
    #ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
    ax.set_yticks(np.arange(0, 700, 200))
    ax.set_yticklabels(np.arange(0, 700, 200), fontsize=14)
    ax.set_ylabel('$N$', fontsize=18)
    ax.text(-12, int(yscale/10)*9.5 + 3, r'$'+txt+'\ \mu_{\\alpha*}$', fontsize=16) 
    ax.text(-12, int(yscale/10)*8.5 + 3, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=16)
    #ax.text(-12, int(yscale/10)*7.5 + 3, "$DEC<6\deg$ ") 
    ax.text(4, int(yscale/10)*9.5 + 3, r"$\bar\mu_{\alpha*}=%0.2f$" %mu, fontsize=16) 
    ax.text(4, int(yscale/10)*8.5 + 3, "$\sigma=%0.2f$" %sigma, fontsize=16)
    ax.set_xlim([xmin, xmax])
    #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
    ax.set_xlabel('')

    Nbin = 15
    ax = axs.ravel()[1]
    histo = ax.hist(muSDSS, Nbin, [xmin, xmax],color='k')
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
    mu2 = np.median(muSDSS)
    #print muSDSS
    sigmaSDSS = (0.741*(np.percentile(muSDSS, 75) - np.percentile(muSDSS, 25)))
    print mu2, sigmaSDSS
    ax.plot([0, 0], [0, 19000], '--', markersize=2.5, color='y')
    ax.set_ylim([0, np.max(histo[0])*1.1])	
    #ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
    ax.set_xticks(np.arange(-10, 12, 5))
    ax.set_yticks(np.arange(0, 700, 200))
    ax.set_yticklabels(np.arange(0, 700, 200), fontsize=14)
    ax.set_xticklabels(np.arange(-10, 12, 5), fontsize=14)
    ax.set_ylabel('$N$', fontsize=18)
    ax.text(-12, int(yscale/10)*9.5, r'$'+txt+'\ \mu_{\delta}$', fontsize=16) 
    #ax.text(-12, int(yscale/10)*7.5, "$DEC<6\deg$ ") 
    ax.text(-12, int(yscale/10)*8.5, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=16) 
    ax.text(4, int(yscale/10)*9.5, r"$\bar\mu_{\delta}=%0.2f$" %mu2, fontsize=16) 
    ax.text(4, int(yscale/10)*8.5, "$\sigma=%0.2f$" %sigmaSDSS, fontsize=16)
    ax.set_xlim([xmin, xmax])
    #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
    ax.set_xlabel(r'$\mu(mas/yr)$', fontsize=18)

    plt.savefig(fname) 

def hist2(x1, x2, txt, root):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(4.5,5.0))
    plt.subplots_adjust(left=0.17, bottom=0.13, right=0.98, top=0.99, wspace=0.16, hspace=0.)
    Nbin = 40
    ax = axs.ravel()[0]  
    histo = ax.hist(x1, Nbin, [0, 4],color='k')
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
    ax.plot([np.median(x1), np.median(x1)], [0, 19000], '--', markersize=2.5, color='y')
    ax.set_ylim([0, np.max(histo[0])*1.1])	
    #ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
    ax.set_ylabel('$N$', fontsize=18)
    ax.text(1.5, int(yscale/10)*8.5,  "$" + txt + "\ Stars$") 
    ax.text(1.5, int(yscale/10)*7.5, "$13.5<r<20.0$")
    ax.text(1.7, int(yscale/10)*4.5, r"$\mu_{\alpha*}$", fontsize=16)
    ax.set_xlim([0,2.6])
    ax.minorticks_on()
    #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
    ax.set_xlabel('')

    ax = axs.ravel()[1]
    histo = ax.hist(x2, Nbin, [0, 4],color='k')
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #ax.plot(histo[1][0:Nbin], histo[0], "r", linewidth=2)
    ax.plot([np.median(x2), np.median(x2)], [0, 19000], '--', markersize=2.5, color='y')
    ax.set_ylim([0, np.max(histo[0])*1.1])	
    #ax.set_yticks([0, int(yscale/10)*4, int(yscale/10)*8])
    ax.set_yticks(np.arange(0, 5000, 1000))
    ax.set_ylabel('$N$', fontsize=18)
    ax.text(1.5, int(yscale/10)*8.5, "$" + txt + "\ Stars$") 
    ax.text(1.5, int(yscale/10)*7.5, "$13.5<r<20.0$") 
    ax.text(1.7, int(yscale/10)*4.5, r"$\mu_{\delta}$", fontsize=16)
    ax.set_xlim([0,2.6])
    ax.minorticks_on()
    #ax.set_yticks([0, 5, 10, 15, 20, 25, 30,35, 40, 45])
    ax.set_xlabel(r'$\chi_{\nu}^{2}$', fontsize=18)
    #ax.set_xlabel(r'$\chi^{2}$', fontsize=18)

    plt.savefig(root+"figs/hist_chi2_"+txt+"_new.pdf") 


def scatterMuCMP(rmag, mu1, mu2, txt, root):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(4.,6.8))
    plt.subplots_adjust(left=0.154, bottom=0.08, right=0.98, top=0.99, wspace=0.16, hspace=0.)
    ax = axs.ravel()[0]
    mask = (rmag>14.0) & (rmag<17.5)

    dltMu = mu1[mask] - mu2[mask]
    print "dltMu.size:", dltMu.size
    avgMuErr = (0.741*(np.percentile(dltMu, 75) - np.percentile(dltMu, 25)))
    print avgMuErr
    subax = plt.axes([0.71,0.60,0.25,0.2])
    histo = subax.hist(dltMu, 20, [-30, 30],color='k')
    #subax.plot(histo[1][0:20], histo[0], "r", linewidth=1.)
    subax.plot([np.median(dltMu),np.median(dltMu)], [0, np.max(histo[0])*1.02], '--', markersize=1.5, color='w')
    subax.set_xlim(-30,30)
    subax.set_xlabel(r'$\Delta\mu_{\alpha}\cos(\delta)(mas/yr)$', fontsize='8')
    subax.set_ylabel(r'$N$', fontsize='8')
    subax.set_xticks([-20, 0, 20])
    subax.set_yticks([0, 200, 400])
    subax.tick_params(labelsize=8)

    plt.axes([0.25,0.88,0.25,0.16])
    minfig = plt.errorbar([0,0], [0,0], avgMuErr, fmt='.', color='r')
    plt.xlim(-0.5,0.5)
    plt.ylim(-100,155)
    plt.axis('off')

    ax.scatter(mu1[mask], mu2[mask],  color='b', marker='.', lw=2.0, alpha=1.0)
    ax.plot([-180, 180], [-180, 180], 'y--', linewidth=1.5)
    #ax.text(-80, 110, "PS1 + SDSS "+txt) 
    ax.text(-80, 90, "14.0$<$r$<$17.5")
    ax.set_xlim([-100,155])
    ax.set_ylim([-100,155])
    ax.set_yticks(np.arange(-90, 136, 45))
    #ax.set_yticks(np.arange(0,1.4,0.25))
    ax.set_xlabel('')
    ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr), OURs$')


    ax = axs.ravel()[1]
    mask = (rmag>17.5) & (rmag<22.5)
    ax.scatter(mu1[mask], mu2[mask],  color='b', marker='.', lw=2.0, alpha=1.0)
    ax.plot([-180, 180], [-180, 180], 'y--', linewidth=1.5)
    #ax.text(-80, 110, "PS1 + SDSS "+txt) 
    ax.text(-80, 90, "17.5$<$r$<$22.5")
    ax.set_xlim([-100,155])
    ax.set_ylim([-100,155])
    ax.set_yticks(np.arange(-90, 136, 45))
    ax.set_xticks(np.arange(-90, 136, 45))
    ax.set_ylabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr), OURs$')
    ax.set_xlabel(r'$\mu_{\alpha}\cos(\delta)(mas/yr), PS1$')

    dltMu = mu1[mask] - mu2[mask]
    print "dltMu.size:", dltMu.size
    avgMuErr = (0.741*(np.percentile(dltMu, 75) - np.percentile(dltMu, 25)))
    print avgMuErr
    subax = plt.axes([0.71,0.14,0.25,0.2])
    histo = subax.hist(dltMu, 20, [-30, 30],color='k')
    #subax.plot(histo[1][0:20], histo[0], "r", linewidth=1.)
    subax.plot([np.median(dltMu),np.median(dltMu)], [0, np.max(histo[0])*1.02], '--', markersize=1.5, color='w')
    subax.set_xlim(-30,30)
    subax.set_xlabel(r'$\Delta\mu_{\alpha}\cos(\delta)(mas/yr)$', fontsize='8')
    subax.set_ylabel(r'$N$', fontsize='8')
    subax.set_xticks([-20, 0, 20])
    subax.set_yticks([0, 200, 400])
    subax.tick_params(labelsize=8)

    plt.axes([0.25,0.425,0.25,0.16])
    minfig = plt.errorbar([0,0], [0,0], avgMuErr, fmt='.', color='r')
    plt.xlim(-0.5,0.5)
    plt.ylim(-100,155)
    plt.axis('off')

    plt.savefig(root+"figs/cmp_ours_PS1_oringinal_"+txt+".png") 

def scatterMuCMP2(rmag, mu1, muErr1, mu2, muErr2, mr0, mr1, txt, flag, root):
    fig = plt.figure(figsize=(5.0,4.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.165, bottom=0.135, right=0.96, top=0.91, wspace=0., hspace=0.)
    mask = (rmag>mr0) & (rmag<mr1)
    if(flag):
    	dltMu = (mu1[mask] - mu2[mask])/np.sqrt(muErr1[mask]**2 + muErr2[mask]**2)
    else:
	dltMu = (mu1[mask] - mu2[mask])
    print "dltMu.size:", dltMu.size
    avgMuErr = (0.741*(np.percentile(dltMu, 75) - np.percentile(dltMu, 25)))
    print np.median(dltMu), avgMuErr
    subax = plt.axes([0.7,0.24,0.235,0.38])
    #subax = plt.axes([0.2,0.21,0.235,0.38])
    histo = subax.hist(dltMu, 30, [-6, 6],color='k')
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #subax.plot(histo[1][0:20], histo[0], "r", linewidth=1.)
    subax.plot([np.median(dltMu),np.median(dltMu)], [0, np.max(histo[0])*1.08], '--', markersize=1.5, color='w')
    subax.set_xlim(-6,6)
    subax.set_ylim([0, np.max(histo[0])*1.17])
    subax.text(-5.6, np.max(histo[0])*1.05, r"$\tilde{\Delta}\bar\mu=%0.2f$" %np.median(dltMu), fontsize='7')
    subax.text(0.9, np.max(histo[0])*1.05, r"$\sigma=%0.2f$" %avgMuErr, fontsize='7')
    subax.set_xlabel(txt[3], fontsize='10') #(mas/yr)
    subax.set_ylabel(r'$N$', fontsize='10')
    subax.set_xticks([-5, 0, 5])
    subax.set_yticks([0, 60, 120]) #([0, int(np.max(histo[0])/3), 2*int(np.max(histo[0])/3)])
    subax.tick_params(labelsize=8)

    '''plt.axes([0.25,0.78,0.25,0.16])
    minfig = plt.errorbar([0,0], [0,0], avgMuErr, fmt='.', color='r')
    plt.xlim(-0.5,0.5)
    plt.ylim(-120,120)
    plt.axis('off')'''

    ax.scatter(mu1[mask], mu2[mask],  color='b', marker='.', lw=2.0, alpha=1.0, rasterized=True)
    ax.plot([-180, 180], [-180, 180], 'y--', linewidth=1.5)
    #ax.text(-80, 110, "PS1 + SDSS "+txt[0]) 
    ax.text(-80, 40, "$"+str(mr0)+"<r<"+str(mr1)+"$", fontsize=14)
    ax.set_xlim([-120,60])
    ax.set_ylim([-120,60])
    ax.set_yticks(np.arange(-100, 60, 50))
    ax.set_xticks(np.arange(-100, 60, 50))
    ax.set_xlabel(txt[1], fontsize=18)
    ax.set_ylabel(txt[2], fontsize=18)#noX


    plt.savefig(root+txt[4]) 

def scatterMuErrCMP(rmag, muErr1, muErr2, txt, root):
    fig = plt.figure(figsize=(5.0,4.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.137, bottom=0.14, right=0.98, top=0.96, wspace=0., hspace=0.)
    mask = (rmag>13.5) & (rmag<20.0)
    dltMuErr = (muErr1[mask] - muErr2[mask])
    print "dltMu.size:", dltMuErr.size
    avgMuErr = (0.741*(np.percentile(dltMuErr, 75) - np.percentile(dltMuErr, 25)))
    print np.median(dltMuErr), avgMuErr
    subax = plt.axes([0.65,0.26,0.3,0.5])
    #subax = plt.axes([0.2,0.21,0.235,0.38])
    histo = subax.hist(dltMuErr, 30, [0, 2.0],color='k') #[0, 4.0]
    yscale = np.min(histo[0])+(np.max(histo[0])-np.min(histo[0]))
    #subax.plot(histo[1][0:20], histo[0], "r", linewidth=1.)
    subax.plot([np.median(dltMuErr),np.median(dltMuErr)], [0, np.max(histo[0])*1.08], '--', markersize=1.5, color='w')
    subax.set_xlim(0, 1.3)#(0,2.3)
    subax.set_ylim([0, np.max(histo[0])*1.17])
    subax.text(0.1, np.max(histo[0])*1.05, r"$\mu=%0.2f$" %np.median(dltMuErr), fontsize='9')#0.2
    subax.text(0.7, np.max(histo[0])*1.05, r"$\sigma=%0.2f$" %avgMuErr, fontsize='9')#1.2
    subax.set_xlabel(r'$\Delta\mu_{\delta,noGAIA} - \Delta\mu_{\delta,GPS1}(mas/yr)$', fontsize='9')
    subax.set_ylabel(r'$N$', fontsize='9')
    #subax.set_xticks(np.arange(0, 2.3, 0.5)) #(0, 1.3, 0.5)
    subax.set_xticks(np.arange(0, 1.3, 0.5)) #
    subax.set_yticks([0, int(np.max(histo[0])/3), 2*int(np.max(histo[0])/3)])
    subax.tick_params(labelsize=8)


    ax.scatter(muErr2[mask], muErr1[mask],  color='b', marker='.', lw=2.0, alpha=1.0)
    #ax.plot([0.6, 4.2], [0.6, 4.2], 'y--', linewidth=1.5) #0.5, 2.2
    ax.plot([0.5, 2.2], [0.5, 2.2], 'y--', linewidth=1.5) #
    #ax.text(1.3, 2.0, "ONLY PS1") 
    ax.text(0.8, 0.8, "13.0$<$r$<$19.0") #1.2, 0.9
    #ax.set_xlim([0.6,3.2]) #0.5,2.1
    #ax.set_ylim([0.6,3.2])
    #ax.set_yticks(np.arange(1.0, 3.2, 0.5)) #0.6, 2.2, 0.2
    #ax.set_xticks(np.arange(1.0, 3.2, 0.5))
    ax.set_xlim([0.6,2.1]) #
    ax.set_ylim([0.6,2.1])
    ax.set_yticks(np.arange(0.8, 2.2, 0.2)) #
    ax.set_xticks(np.arange(0.8, 2.2, 0.2))
    ax.set_ylabel(r'$\Delta\mu_{\delta}(mas/yr), no GAIA$')#PS1$')
    ax.set_xlabel(r'$\Delta\mu_{\delta}(mas/yr), GPS1$')
    #plt.show()
    plt.savefig(root+"cmp_muErr_GPS1_300.eps") 

def cmpObjs(ra0, dec0, ra, dec, root):
    uniqueStar = np.unique(obj_id)
    print "len(uniqueStar)", len(uniqueStar)
    R2 = []
    idx=0
    count = 0
    for idx in range(0, 1):#len(uniqueStar)):
	ind = obj_id == -6777811886076048565#uniqueStar[idx]
	xtmp = mjd[ind]
	tmp = ra[ind]
	ytmp = (tmp-np.min(tmp))*1000*3600
	yErrtmp = raErr[ind]#*1000*3600
	objIDtmp = obj_id[ind][0]
	rMagtmp = mr[ind][0]
	if((len(xtmp)>2) and (np.max(ytmp)>2)):
		fig = plt.figure(figsize=(6.0,4.0))
		ax = fig.add_subplot(111)
    		plt.subplots_adjust(left=0.1, bottom=0.11, right=0.96, top=0.91, wspace=0., hspace=0.)
		ax.errorbar(xtmp, ytmp, yErrtmp, fmt='.', color='red')

def plotXvalidation(dxX, ErrDec, nObs, root):
    plt.figure()
    plt.hist(dxX, 40, [-45, 45])
    plt.xlabel(r'$\Delta \delta (mas)$')
    plt.ylabel('N')
    plt.xlim([-6,6.0])
    plt.text(-5, 1400, r"$\langle\Delta\delta/Err_{\delta}\rangle=%0.2f$" %np.median(dxX))
    plt.text(-5, 1340, r"$\sigma=%0.2f$" %(0.741*(np.percentile(dxX, 75) - np.percentile(dxX, 25))))
    #plt.show()
    plt.savefig(root+"figs/hist_Ddec_ratio.png")

    plt.figure()
    plt.hist(ErrDec, 40, [-30, 30])
    plt.xlabel(r'$Err_{\delta}(mas)$')
    plt.ylabel('N')
    plt.xlim([1,20.0])
    plt.text(15, 2500, r"$\langle Err_{\delta}\rangle=%0.2f$" %np.median(ErrDec))
    plt.text(15, 2350, r"$\sigma=%0.2f$" %(0.741*(np.percentile(ErrDec, 75) - np.percentile(ErrDec, 25))))
    #plt.show()
    plt.savefig(root+"figs/hist_Errdec.png")

    plt.figure()
    plt.scatter(nObs, dxX, color='b', marker='.', lw=2.0, alpha=0.2)
    binned = med_bin(nObs, dxX, 1.5, 20.5, 1)
    good = binned[2] > 0
    plt.errorbar(binned[0][good], binned[1][good], yerr=binned[2][good], fmt='o', color='r', ecolor='r',elinewidth=2.0)

    plt.ylabel(r'$\Delta \delta/Err_{\delta}(mas)$')
    plt.xlabel('Nobs')
    plt.xlim([0,22])
    plt.ylim([-10,10.0])
    #plt.text(-40, 750, r"$\bar{\Delta\delta}=%0.2f$" %np.median(dxX))
    #plt.text(-40, 700, r"$\sigma=%0.2f$" %(0.741*(np.percentile(dxX, 75) - np.percentile(dxX, 25))))
    #plt.show()
    plt.savefig(root+"figs/Nobs_Ddec.png")

def plotDDY(x_fix, y_fix, ddy, root):

    nbin_x = 32.
    nbin_y = 32.
    x_min = 0.
    x_max = 4096.
    y_min = 0.
    y_max = 4096.
    good = np.where((x_fix >= x_min) & (x_fix <= x_max) & (y_fix >= y_min) & (y_fix <= y_max))[0]
    h = esutil.stat.histogram2d(x_fix[good], y_fix[good], nx=nbin_x, ny=nbin_y, xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, rev=True, more=True)

    # calculate the clipped average of each pixel
    med_of_pixel, rms_of_pixel, n_per_pixel = useful.stats_on_pixels2(h, ddy[good], (np.zeros(len(ddy))+1.0)[good])

    # reshape into 2D arrays
    med_of_pixel = med_of_pixel.reshape(nbin_y, nbin_x)
    rms_of_pixel = rms_of_pixel.reshape(nbin_y, nbin_x)
    n_per_pixel = n_per_pixel.reshape(nbin_y, nbin_x)

    # normalize the ZPVM to the central pixel
    med_of_pixel -= med_of_pixel[16,16]

    # find empty pixels
    bad = np.where(med_of_pixel < -90)
    med_of_pixel[bad] = (med_of_pixel[bad[0]+1, bad[1]]+med_of_pixel[bad[0]-1, bad[1]])/2.
    if bad[0].size > 0:
        print 'Interpolated %d empty pixel(s).' % (bad[0].size)


    # if a ZPVM already exist, add the new residuals to it
    '''try:
        with open(root+'zpvm_%s_%s.npy' % (sys.argv[1], sys.argv[2])):
            zpvm = np.load(root+'zpvm_%s_%s.npy' % (sys.argv[1], sys.argv[2]))
            zpvm += med_of_pixel
            np.save('zpvm_%s_%s.npy' % (sys.argv[1], sys.argv[2]), zpvm)
    except IOError:
        # if not, just dump the ZPVM
        np.save(root+'zpvm_%s_%s.npy' % (sys.argv[1], sys.argv[2]), med_of_pixel)'''

    # plot
    good = np.where(med_of_pixel > -90)
    range = np.max([np.abs(med_of_pixel[good].max()), np.abs(med_of_pixel[good].min())])
    print range
    a = 3.47*2
    fig1 = plt.figure(figsize=(a,a), num=1)
    ax = fig1.add_subplot(111, aspect='equal')
    plt.clf()
    cax = plt.imshow(ma.masked_where(med_of_pixel < -90, med_of_pixel), extent=[h['xmin'], h['xmax'], h['ymin'], \
	h['ymax']], origin='lower', cmap=cm.jet, interpolation='nearest', vmin=-1.*range, vmax=range)
    plt.xlabel('x')
    plt.ylabel('y')
    cbar = plt.colorbar(cax, orientation='horizontal')
    #cbar.ax.set_xticklabels(['-0.008','-0.004', '0', '0.004','0.008'], fontsize='medium')
    cbar.ax.set_xlabel('residuals', fontsize='medium')
    plt.savefig(root+'figs/res_xy_000.png')


    # plot
    good = np.where(rms_of_pixel > -90)
    range = rms_of_pixel[good].max()
    print range/2.
    a = 3.47*2
    fig1 = plt.figure(figsize=(a,a), num=1)
    ax = fig1.add_subplot(111, aspect='equal')
    plt.clf()
    cax = plt.imshow(ma.masked_where(rms_of_pixel < -90, rms_of_pixel), extent=[h['xmin'], h['xmax'], h['ymin'], h['ymax']], \
	origin='lower', cmap=cm.jet, interpolation='nearest', vmin=0, vmax=range)
    plt.xlabel('x')
    plt.ylabel('y')
    cbar = plt.colorbar(cax, orientation='horizontal')
    #cbar.ax.set_xticklabels(['-0.008','-0.004', '0', '0.004','0.008'], fontsize='medium')
    cbar.ax.set_xlabel('rms scatter', fontsize='medium')
    plt.savefig(root+'figs/rms_xy_000.png')

    # plot
    a = 3.47*2
    fig2 = plt.figure(figsize=(a,a), num=2)
    ax = fig2.add_subplot(111, aspect='equal')
    plt.clf()
    cax = plt.imshow(n_per_pixel, extent=[h['xmin'], h['xmax'], h['ymin'], h['ymax']], origin='lower', cmap=cm.jet, interpolation='nearest')
    plt.xlabel('x')
    plt.ylabel('y')
    step = int(np.floor((n_per_pixel.max()-n_per_pixel.min())/5.))
    cbar = plt.colorbar(cax, ticks=np.arange(n_per_pixel.min(), n_per_pixel.max(), step), orientation='horizontal')
    cbar.ax.set_xlabel('Sources per pixel', fontsize='medium')
    plt.savefig(root+'figs/counts_xy_000.png')


    chi = ma.masked_where(med_of_pixel < -90, med_of_pixel).flatten()/ma.masked_where(rms_of_pixel < -90, rms_of_pixel).flatten()
    plt.clf()
    plt.hist(chi, 30, [-0.6, 0.6])
    plt.text(-0.5, 100, r"$\mu=%0.2f$" %np.median(chi)) #\langle RES/RMS\rangle
    plt.text(-0.5, 95, r"$\sigma=%0.2f$" %(0.741*(np.percentile(chi, 75) - np.percentile(chi, 25))))
    plt.xlabel(r'$RES/RMS$')
    plt.ylabel(r'$N$')
    plt.xlim([-0.6,0.6])
    plt.ylim([0,110])
    plt.savefig(root+'figs/hist_RES_RMS_000.png')

    '''fig=plt.figure()
    #cax1 = fig.add_axes([0.3963, 0.90, 0.2955, 0.016])
    cmap1 = mpl.cm.jet  #cool
    norm1 = mpl.colors.Normalize(vmin=-50, vmax=50)#(vmin=0, vmax=7.5)
    plt.scatter(x_fix, y_fix, c=ddy, norm = norm1, cmap=cmap1, marker='.', lw=0)
    plt.xlabel('$X_{FIX}$')
    plt.ylabel('$Y_{FIX}$')
    plt.xlim([0,4800])
    plt.ylim([0,4800])
    plt.colorbar()
    #plt.show()
    plt.savefig(root+"figs/offset_Pal5.png")'''
    
##############plot DetaRA v.s. epoch(MJD) for different patches in one panel################
#def xxx():
    '''plt.figure()    
    for pidx in range(0, pN):
    
        #plt.plot(mjdValues, residualRaArray, '.', markersize=1.0, color='k')
        plt.errorbar(medianMJD[pidx, :], medianRAresidual[pidx, :], stdRAresidual[pidx, :], \
            markersize=5.5, color=clor[0], fmt='--o', ecolor=clor[pidx],   capthick=2,elinewidth=2)
        if(pidx==0): plt.plot([0, 60000], [0, 0], '--', markersize=2.5, color='g') 
    plt.xlabel('MJD')
    plt.ylabel('$\Delta RA$/mas')
    plt.xlim(55000,56600)
    plt.ylim(-8.0,7.)
    #plt.xticks(np.array([8, 10, 12, 14]), fontsize='small')
    plt.yticks(np.array([-6, -3, 0, 3, 6]), fontsize='small')
    plt.title('Residual distribution in different Epochs') 
    #plt.grid(True)
    plt.savefig("residualRA"+str(999)+".png")'''



def sm_hist2(data, delta=5):
	dataMin = np.floor(data.min())
	dataMax = np.ceil(data.max())
	n_bin = np.ceil(1.*(dataMax-dataMin) / delta) + 1
	idxs = ((data  - dataMin) / delta).astype(int)
	counts = np.zeros(n_bin) 
	bin_edges = np.arange(dataMin, dataMax+delta, delta)
	for idx in idxs:
		counts[idx] += 1
	#print counts
	# These two lines double the points let you make the histogram
	counts = np.ravel(zip(counts,counts)) 
	bin_edges = np.ravel(zip(bin_edges,bin_edges))
	counts = np.hstack((np.array([0]), counts))
	bin_edges = np.hstack((bin_edges, bin_edges[-1]))
	return counts, bin_edges


		
def plotTable(stars):	
    latex_txt = "\\"
    for idx2 in range(0, len(stars)):		
	latex_txt = str(idx2+1) + "&" + (stars['name'][idx2]) + "&" + str('%.3f'%(stars['ra'][idx2])) + "&" + str('%.3f'%(stars['dec'][idx2]))+"&"
	latex_txt = latex_txt + str('%.2f'%stars['mura'][idx2]) + "$\pm$" + str('%.2f'%stars['muraerr'][idx2]) + "&" \
	    + str('%.2f'%stars['mudec'][idx2]) + "$\pm$" + str('%.2f'%stars['mudecerr'][idx2]) + "&" \
	    + str('%.2f'%stars['gmag'][idx2]) \
	    + "\\" + "\\"			
	print(latex_txt)

def plotTable1(pm):	
    latex_txt = "\\"
    for idx2 in range(0, len(pm)):		
	latex_txt = str(pm['obj_id'][idx2]) + "&" + str('%.3f'%(pm['ra'][idx2])) + "&" + str('%.3f'%(pm['dec'][idx2]))+"&"
	latex_txt = latex_txt + str('%.2f'%pm['muXI'][idx2]) + "$\pm$" + str('%.2f'%pm['muXIerr'][idx2]) + "&" \
	    + str('%.2f'%pm['muETA'][idx2]) + "$\pm$" + str('%.2f'%pm['muETAerr'][idx2]) + "&" \
	    + str('%.2f'%pm['mr'][idx2]) + "&" + str('%.2f'%pm['chi2XI'][idx2]) + "&" + str('%.2f'%pm['chi2ETA'][idx2]) + "&" \
	    +  str(pm['flag'][idx2]) + "\\" + "\\"			
	print(latex_txt)
def plotTable2(pm):	
    latex_txt = "\\"
    for idx2 in range(0, len(pm)):		
	latex_txt = str(pm['objID'][idx2]) + "&" + str('%.3f'%(pm['ra'][idx2])) + "&" + str('%.3f'%(pm['dec'][idx2]))+"&"
	latex_txt = latex_txt + str('%.2f'%pm['muXxi'][idx2]) + "$\pm$" + str('%.2f'%pm['muErrXI'][idx2]) + "&" \
	    + str('%.2f'%pm['muXeta'][idx2]) + "$\pm$" + str('%.2f'%pm['muErrETA'][idx2]) + "&" \
	    + str('%.2f'%pm['mr'][idx2]) + "&" + str('%.2f'%pm['chi2XI'][idx2]) + "&" + str('%.2f'%pm['chi2ETA'][idx2]) + "&" \
	    +  str(pm['Nobs'][idx2]) + "&" +  str(pm['flag'][idx2]) + "\\" + "\\"			
	print(latex_txt)
