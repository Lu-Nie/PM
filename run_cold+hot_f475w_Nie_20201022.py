import os, subprocess, shlex, datetime, multiprocessing as mp, numpy as np, astropy.io.fits as pyfits
import pandas as pd

start_t = datetime.datetime.now()
num_cores = int(mp.cpu_count())
#print("The computer has " + str(num_cores) + " cores")
 
data_path, res_path = './data/', './res_test/'  

#cutting
#os.system('cd '+data_path+'; astcrop --section=1691:1865,698:760 --mode=img --hdu=0 ds9.fits') 
#pref = 'hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f475w_v1_'#f160w
    
#detimage = data_path+pref+'drz_cropped.fits'
#image = data_path+pref+'drz_cropped.fits'
#detweight = data_path+pref+'rms.fits'
#weight = data_path+pref+'rms.fits'

#pref = 'dsB09F01-3.0_'#f475w
 
#detimage = data_path+pref+'drz.fits'
#image = data_path+pref+'drz.fits'
#detweight = data_path+pref+'rms.fits'
#weight = data_path+pref+'rms.fits'

#os.system('cd '+data_path+'; ds9 -width 1360 -height 768  hlsp_phat_hst_acs-wfc_12057-m31-b09-f01_f475w_v1_drz.fits -zoom 4 -pan to 161.875 58.75 image -saveimage ds9b09f01f475w1.fits -exit') 

pref = 'ds9b09f01f814w0'#f160w

detimage = data_path+pref+'.fits'
image = data_path+pref+'.fits'

cold, hot = 'ACS-WFC.Hfinal.rms.cold.Nie2020.sex', 'ACS-WFC.Hfinal.rms.rome.Nie2020.sex'

coldcat, hotcat = res_path+pref+'cold.cat', res_path+pref+'hot.cat'
coldseg, hotseg = res_path+pref+'segcold.fits', res_path+pref+'seghot.fits'
coldaper, hotaper = res_path+pref+'apercold.fits', res_path+pref+'aperhot.fits'

#outcat = res_path+pref+'comb.cat'
outcat = res_path+pref+'comb.cat'
out_bad_cat = res_path+pref+'_badcomb.cat'
outseg = res_path+pref+'segcomb.fits'
outparam = 'default.param'

gain = '2.0'   # [ 7200, 5450, 7028, 18232, 5017.606, 8197, 5086, 5250 ] for
                # [ B,    V,    i,    z,     F098M,    Y,    J,    H    ]
                #[2.0, 1.5] for [ACS, UVIS]
magzp = '25.94' # [ 25.673, 26.486, 25.654, 24.862, 25.68, 26.27, 26.25, 25.96 ] for
                # [ B,      V,      i,      z,      F098M, Y,     J,     H     ]
                #[25.94, 26.05] for [f814w, f475w]
        
seeing = '0.05'#[0.05, 0.09] for [f814w, f475w]
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
if len(cold_table['KRON_RADIUS'][idx_c]) > 0:
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


#####fits image
#image_hdu = []
# for j in range(1,4):
#     seghd, segim = hc[j].header, hc[j].data
#     seghd_hot, segim_hot = hh[j].header, hh[j].data
#     for i in range(0, nhot):
#         print(nhot, i)
#         idx = np.where(segim_hot == hot_table['NUMBER'][i])
#     #print(len(idx), idx, segim[idx])
#         if len(idx) > 1:
#             segim[idx] = off + i
#     #print(segim[idx])
#         fs.write(str(int(off+i))+' '+str(int(hot_table['NUMBER'][i]))+'\n')
#         hot_table['NUMBER'][i] = off + i
#         image_hdu.append(pyfits.ImageHDU(segim))
# fs.close()
# #print(segim)

# primary_hdu = pyfits.PrimaryHDU(header=seghd)
# #image_hdu = pyfits.ImageHDU(segim)
# hdul = pyfits.HDUList([primary_hdu, image_hdu[0], image_hdu[1], image_hdu[2]])
# hdul.writeto(outseg, overwrite=True)


######end
table_all = np.append(cold_table, hot_table)
np.savetxt(outcat, table_all, fmt="%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %d %d %d %d %f %f %d", header='NUMBER FLUX_ISO FLUXERR_ISO MAG_ISO MAGERR_ISO FLUX_AUTO FLUXERR_AUTO MAG_AUTO MAGERR_AUTO MAG_BEST MAGERR_BEST KRON_RADIUS BACKGROUND X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 CXX_IMAGE CYY_IMAGE CXY_IMAGE ELLIPTICITY FWHM_IMAGE FLUX_RADIUS FLAGS CLASS_STAR A_IMAGE B_IMAGE XMIN_IMAGE YMIN_IMAGE XMAX_IMAGE YMAX_IMAGE ELONGATION THETA_IMAGE ISOAREA_IMAGE')

reg1 = res_path + 'test_f475w1.0.reg'
reg2 = res_path + 'test_f475w2.0.reg'
os.system("awk '{print \"ellipse(\"$14\",\"$15\",\"($26*$12)\",\"($27*$12)\",\"$33\")\"}' "+outcat+ " > " + reg1)
os.system("awk '{print $16,$17}' "+outcat+ " > " + reg2)

#psf


end_t = datetime.datetime.now()
elapsed_sec = (end_t - start_t).total_seconds()
print("Used Time: " + "{:.2f}".format(elapsed_sec) + " sec")
