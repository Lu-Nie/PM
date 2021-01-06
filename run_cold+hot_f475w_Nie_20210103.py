import os, subprocess, shlex, datetime, multiprocessing as mp, numpy as np, astropy.io.fits as pyfits
import pandas as pd

start_t = datetime.datetime.now()
num_cores = int(mp.cpu_count())
#print("The computer has " + str(num_cores) + " cores")
 
data_path, res_path = './data/', './res_test/'  


X = np.arange(323.75,3824.25,647.5)
Y = np.arange(117.5,1191.5,235)
N = 0
table_all2 = []
outcat = res_path+'ds9b09f01f814w_comb.cat'
for p in range(len(Y)):
    for q in range(len(X)):
        x = X[q]
        y = Y[p]
        n = str(N)
        os.system('cd '+data_path+'; ds9 -width 1360 -height 768 hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f814w_v1_drz.fits -zoom 2 -pan to %f %f image -saveimage ds9b09f01f814w%s.fits -exit' %(x,y,n))
        
        pref = 'ds9b09f01f814w'+ n #f475w

        detimage = data_path+pref+'.fits'
        image = data_path+pref+'.fits'

        cold, hot = 'ACS-WFC.Hfinal.rms.cold.Nie2020.sex', 'ACS-WFC.Hfinal.rms.rome.Nie2020.sex'

        coldcat, hotcat = res_path+pref+'cold'+n+'.cat', res_path+pref+'hot'+n+'.cat'
        coldseg, hotseg = res_path+pref+'segcold.fits', res_path+pref+'seghot.fits'
        coldaper, hotaper = res_path+pref+'apercold.fits', res_path+pref+'aperhot.fits'

        #outcat = res_path+pref+'comb.cat'
        #outcat = res_path+pref+'comb.cat'
        out_bad_cat = res_path+pref+'_badcomb.cat'
        outseg = res_path+pref+'segcomb.fits'
        outparam = 'default.param'

        gain = '2.0'   # [ 7200, 5450, 7028, 18232, 5017.606, 8197, 5086, 5250 ] for
                        # [ B,    V,    i,    z,     F098M,    Y,    J,    H    ]
        magzp = '26.05' # [ 25.673, 26.486, 25.654, 24.862, 25.68, 26.27, 26.25, 25.96 ] for
                        # [ B,      V,      i,      z,      F098M, Y,     J,     H     ]
        #seeing = '0.18'
        seeing = '0.09'
        #DETth_hot = '0.3'
        hotidtable = res_path+pref+'hotid.cat'

        #--------------------------------------
        # Run cold 
        print('Running cold')
        os.system("sex "+detimage+","+image+" -c "+cold+" -CATALOG_NAME "+coldcat+" -CATALOG_TYPE ASCII "+\
        " -PARAMETERS_NAME "+outparam+" -WEIGHT_TYPE NONE,NONE -CHECKIMAGE_TYPE SEGMENTATION,APERTURES -CHECKIMAGE_NAME "+\
        coldseg+","+coldaper+" -GAIN "+gain+" -MAG_ZEROPOINT "+magzp+" -SEEING_FWHM "+seeing)

        # Run hot
        print('Running hot')
        os.system("sex "+detimage+","+image+" -c "+hot+" -CATALOG_NAME "+hotcat+" -CATALOG_TYPE ASCII "+\
        " -PARAMETERS_NAME "+outparam+" -WEIGHT_TYPE NONE,NONE -CHECKIMAGE_TYPE SEGMENTATION,APERTURES -CHECKIMAGE_NAME "+\
        hotseg+","+hotaper+" -GAIN "+gain+" -MAG_ZEROPOINT "+magzp+" -SEEING_FWHM "+seeing)
        #+" -DETECT_THRESH "+DETth_hot
        """
        # no check images
        # Run cold 
        print('Running cold')
        os.system("sex "+detimage+","+image+" -c "+cold+" -CATALOG_NAME "+coldcat+" -CATALOG_TYPE ASCII "+\
        " -PARAMETERS_NAME "+outparam+" -WEIGHT_IMAGE "+detweight+","+weight+" -WEIGHT_TYPE MAP_RMS,MAP_RMS -CHECKIMAGE_TYPE NONE "+\
        " -GAIN "+gain+" -MAG_ZEROPOINT "+magzp+" -SEEING_FWHM "+seeing)

        # Run hot
        print('Running hot')
        os.system("sex "+detimage+","+image+" -c "+hot+" -CATALOG_NAME "+hotcat+" -CATALOG_TYPE ASCII "+\
        " -PARAMETERS_NAME "+outparam+" -WEIGHT_IMAGE "+detweight+","+weight+" -WEIGHT_TYPE MAP_RMS,MAP_RMS -CHECKIMAGE_TYPE NONE "+\
        " -GAIN "+gain+" -MAG_ZEROPOINT "+magzp+" -SEEING_FWHM "+seeing)
        """

        #--------------------------------------
        # Read hotcat and coldcat 
        print('Read cold and hot catalogs')
        a = open(outparam,'r').read().split('\n')
        h = [item for item in a if item!='' and item[0]!='#']
        print(h)

        cold_table = np.genfromtxt(coldcat, names=h) # 22223
        idx_c = np.where(cold_table['KRON_RADIUS']==0)
        if len(idx_c) > 1 and len(cold_table['KRON_RADIUS'][idx_c]) > 0:
            cold_table['KRON_RADIUS'][idx_c] = np.median(cold_table['KRON_RADIUS'])
        #print(len(cold_table['KRON_RADIUS'][idx_c])) # 0

        hot_table = np.genfromtxt(hotcat, names=h) # 39428
        idx_h = np.where(hot_table['KRON_RADIUS']==0)
        if len(hot_table['KRON_RADIUS'][idx_h]) > 0:
            hot_table['KRON_RADIUS'][idx_h] = np.median(hot_table['KRON_RADIUS'])
        #print(len(hot_table['KRON_RADIUS'][idx_h])) # 62

        #--------------------------------------
        print('Including hot detections')
        cold_cxx = cold_table['CXX_IMAGE'] / cold_table['KRON_RADIUS']**2 # 22223
        cold_cyy = cold_table['CYY_IMAGE'] / cold_table['KRON_RADIUS']**2
        cold_cxy = cold_table['CXY_IMAGE'] / cold_table['KRON_RADIUS']**2
        ncold = len(cold_table) # ncold = 22223

        hot_cxx = hot_table['CXX_IMAGE'] / hot_table['KRON_RADIUS']**2
        hot_cyy = hot_table['CYY_IMAGE'] / hot_table['KRON_RADIUS']**2
        hot_cxy = hot_table['CXY_IMAGE'] / hot_table['KRON_RADIUS']**2
        nhot = len(hot_table) # nhot = 39428

        hc = pyfits.open(coldseg)
        seghd, segim = hc[0].header, hc[0].data# [40500, 32400]

        hh = pyfits.open(hotseg)
        seghd_hot, segim_hot = hh[0].header, hh[0].data # [40500, 32400]

        #------------------------------------------

        for i in range(0, ncold): # range(0, ncold) # ncold = 22223
            print('N:', ncold, i, len(cold_cxx[i] * (hot_table['X_IMAGE'] - cold_table['X_IMAGE'][i])**2 + cold_cyy[i] * (hot_table['Y_IMAGE'] - cold_table['Y_IMAGE'][i])**2 + cold_cxy[i] * (hot_table['X_IMAGE'] - cold_table['X_IMAGE'][i]) * (hot_table['Y_IMAGE'] - cold_table['Y_IMAGE'][i])))
            idx = np.where( cold_cxx[i] * (hot_table['X_IMAGE'] - cold_table['X_IMAGE'][i])**2 + \
                            cold_cyy[i] * (hot_table['Y_IMAGE'] - cold_table['Y_IMAGE'][i])**2 + \
                            cold_cxy[i] * (hot_table['X_IMAGE'] - cold_table['X_IMAGE'][i]) * (hot_table['Y_IMAGE'] - cold_table['Y_IMAGE'][i]) > 1.1**2 )
            hot_table = hot_table[idx]
        print(len(hot_table)) # 14091

        #--------------------------------------
        # Read the segmentaiton images and add objects from hot segmentation map to cold segementation map,
        # but only at pixels where no object was defined in cold segmentation map. Then write result.
        print('Creating combined segmentation map')

        os.system('> '+hotidtable)
        fs = open(hotidtable, 'a')

        nhot = len(hot_table) # 14091
        off = np.max(cold_table['NUMBER']) + 1

        for i in range(0, nhot):
            print(nhot, i)
            idx = np.where(segim_hot == hot_table['NUMBER'][i])
            #print(len(idx), idx, segim[idx])
            if len(idx) > 0:
                segim[idx] = off + i
            #print(segim[idx])
            fs.write(str(int(off+i))+' '+str(int(hot_table['NUMBER'][i]))+'\n')
            hot_table['NUMBER'][i] = off + i
        fs.close()

        primary_hdu = pyfits.PrimaryHDU(header=seghd)
        image_hdu = pyfits.ImageHDU(segim)
        hdul = pyfits.HDUList([primary_hdu, image_hdu])
        hdul.writeto(outseg, overwrite=True)
        N = N+1
        table_all = np.append(cold_table,hot_table)
        table_all2.append(table_all)
        if N > 0:
            break
table_all3 = table_all2[0]
for j in range(1,len(table_all2)):
    table_all3  = np.append(table_all3,table_all2[j])
np.savetxt(outcat, table_all3, fmt="%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %d %d %d %d %f %f %d", header='NUMBER FLUX_ISO FLUXERR_ISO MAG_ISO MAGERR_ISO FLUX_AUTO FLUXERR_AUTO MAG_AUTO MAGERR_AUTO MAG_BEST MAGERR_BEST KRON_RADIUS BACKGROUND X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 CXX_IMAGE CYY_IMAGE CXY_IMAGE ELLIPTICITY FWHM_IMAGE FLUX_RADIUS FLAGS CLASS_STAR A_IMAGE B_IMAGE XMIN_IMAGE YMIN_IMAGE XMAX_IMAGE YMAX_IMAGE ELONGATION THETA_IMAGE ISOAREA_IMAGE')

reg1 = res_path + 'test_f475w1.0.reg'
reg2 = res_path + 'test_f475w2.0.reg'
os.system("awk '{print \"ellipse(\"$14\",\"$15\",\"($26*$12)\",\"($27*$12)\",\"$33\")\"}' "+outcat+ " > " + reg1)
os.system("awk '{print $16,$17}' "+outcat+ " > " + reg2)

#psf


end_t = datetime.datetime.now()
elapsed_sec = (end_t - start_t).total_seconds()
print("Used Time: " + "{:.2f}".format(elapsed_sec) + " sec")
