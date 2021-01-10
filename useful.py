import numpy as np
import scipy.interpolate
import matplotlib
import math
#import baryvel
from scipy.integrate import simps
import sys
from esutil.coords import eq2gal
from numpy.linalg import inv

#extinction_coeff = [1.87, 1.38, 1.00, 0.76, 0.54, 0.327, 0.209, 0.133] # ugrizJHK
extinction_coeff = [1.810, 1.400, 1.000, 0.759, 0.561, 0.317, 0.200, 0.132] # ugrizJHK according to Berry et al. (2012)

def degToHMS(deg):
    hrs = np.int(np.floor(deg/15.0))
    temp = (deg/15.0 - hrs)*60.0
    min = np.int(np.floor(temp))
    sec = (temp - min)*60.0
    return (hrs, min, sec)

def degToDMS(deg):
    s = 1.0
    if(deg < 0.0):
        s = -1.0
        deg *= -1.0
    degrees = int(math.floor(deg))
    temp = (deg - degrees)*60.0
    min = int(math.floor(temp))
    sec = (temp - min)*60.0
    return ("+" if s >= 0.0 else "-", abs(degrees), abs(min), abs(sec))

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = np.linspace(0,1.,N)
    # N+1 indices
    indices = np.linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = np.array(cdict[key])
        I = scipy.interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = np.zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def stats_on_pixels(h, dat, min_number_of_points=5):
    """
    Description: Find median for each pixel
    Example: med_of_pixel, count_in_pixel = stats_on_pixels(h, master['stdev'])
    h = dictionary returned by esutil's histogram2d with more=True
    dat = array on which the statistics will be calculated
    """
    # loop over all bins
    n_pixels = h['rev'][0]-1
    pixel_list = list()
    for j in np.arange(n_pixels):
        if (h['rev'][j+1] > h['rev'][j]):
            pixel_array = h['rev'][h['rev'][j]:h['rev'][j+1]]
        else:
            pixel_array = np.array([-1])
        pixel_list.append(pixel_array)

    count_in_pixel = []
    med_of_pixel = []
    for j in np.arange(n_pixels):
        pixel = pixel_list[j]
        if pixel[0] != -1:
            count_in_pixel = np.append(count_in_pixel, np.int(pixel.size))
            if pixel.size >= min_number_of_points:
                med_of_pixel = np.append(med_of_pixel, np.median(dat[pixel]))
            else:
                med_of_pixel = np.append(med_of_pixel, -99.99)
        else:
            count_in_pixel = np.append(count_in_pixel, 0)
            med_of_pixel = np.append(med_of_pixel, -99.99)
    return med_of_pixel, count_in_pixel

def sum_in_pixels(h, dat, min_number_of_points=1):
    """
    Description: Find the sum in each pixel
    Example: sum_in_pixel, count_in_pixel = sum_in_pixels(h, LF_all)
    h = dictionary returned by esutil's histogram2d with more=True
    dat = array on which the statistics will be calculated
    """
    # loop over all bins
    n_pixels = h['rev'][0]-1
    pixel_list = list()
    for j in np.arange(n_pixels):
        if (h['rev'][j+1] > h['rev'][j]):
            pixel_array = h['rev'][h['rev'][j]:h['rev'][j+1]]
        else:
            pixel_array = np.array([-1])
        pixel_list.append(pixel_array)

    count_in_pixel = []
    sum_in_pixel = []
    for j in np.arange(n_pixels):
        pixel = pixel_list[j]
        if pixel[0] != -1:
            count_in_pixel = np.append(count_in_pixel, np.int(pixel.size))
            if pixel.size >= min_number_of_points:
                sum_in_pixel = np.append(sum_in_pixel, np.sum(dat[pixel]))
            else:
                sum_in_pixel = np.append(sum_in_pixel, 0.)
        else:
            count_in_pixel = np.append(count_in_pixel, 0)
            sum_in_pixel = np.append(sum_in_pixel, 0.)
    return sum_in_pixel, count_in_pixel

def stats_on_pixels2(h, dat, datErr, min_number_of_points=5):
    """
    Description: Find clipped, weighted average for each pixel
    Example: avg_of_pixel, count_in_pixel = stats_on_pixels2(h, dat['res'], dat['resErr'])
    h = dictionary returned by esutil's histogram2d with more=True
    dat = array on which the statistics will be calculated
    """
    # loop over all bins
    n_pixels = h['rev'][0]-1
    pixel_list = list()
    for j in np.arange(n_pixels):
        if (h['rev'][j+1] > h['rev'][j]):
            pixel_array = h['rev'][h['rev'][j]:h['rev'][j+1]]
        else:
            pixel_array = np.array([-1])
        pixel_list.append(pixel_array)

    count_in_pixel = []
    avg_of_pixel = []
    rms_of_pixel = []
    for j in np.arange(n_pixels):
        pixel = pixel_list[j]
        if pixel[0] != -1:
            if pixel.size >= min_number_of_points:
                rms = 0.741*(np.percentile(dat[pixel], 75) - np.percentile(dat[pixel], 25))
                good = np.where(np.abs(np.median(dat[pixel])-dat[pixel]) < 3*rms)[0]
                if good.size >= min_number_of_points:
                    avg_of_pixel = np.append(avg_of_pixel, np.average(dat[pixel[good]], weights=1./datErr[pixel[good]]**2))
                    rms_of_pixel = np.append(rms_of_pixel, np.std(dat[pixel[good]], ddof=1))
                else:
                    avg_of_pixel = np.append(avg_of_pixel, -99.99)
                    rms_of_pixel = np.append(rms_of_pixel, -99.99)
                count_in_pixel = np.append(count_in_pixel, np.int(good.size))
            else:
                avg_of_pixel = np.append(avg_of_pixel, -99.99)
                rms_of_pixel = np.append(rms_of_pixel, -99.99)
                count_in_pixel = np.append(count_in_pixel, np.int(good.size))
        else:
            avg_of_pixel = np.append(avg_of_pixel, -99.99)
            rms_of_pixel = np.append(rms_of_pixel, -99.99)
            count_in_pixel = np.append(count_in_pixel, 0)
    return avg_of_pixel, rms_of_pixel, count_in_pixel

def rebin(a, *args):
    '''
    rebin ndarray data into a smaller ndarray of the same rank whose
    dimensions are factors of the original dimensions. eg. An array with 6
    columns and 4 rows can be reduced to have 6,3,2 or 1 columns and 4,2 or 1
    rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

def EarthBaryVel(jd):
    """
    Description: Returns barycentric velocity components of Earth given a Julian date
    Example: bvx, bvy, bvz = EarthBaryVel(jd)
    """
    result = baryvel.baryvel(jd)
    bvx, bvy, bvz = result[1]
    return bvx, bvy, bvz

def bary_vel(rv, ra, dec, jd):
    """
    Description: Provides barycentric velocity.
    Example: bary_rv = bary_vel(rv, ra, dec, jd)
    """
    bvx, bvy, bvz = EarthBaryVel(jd)
    ra = np.radians(ra)
    dec = np.radians(dec)
    return rv+bvx*np.cos(dec)*np.cos(ra) + bvy*np.cos(dec)*np.sin(ra) + bvz*np.sin(dec)

baryVel = np.vectorize(bary_vel)


def gal2XYZ(gl, gb, helio_dist, Rsun=8.00):
    """
    Description: Convert J2000.0 equatorial RA, Dec, and heliocentric distance in kpc to Galactocentric coordinates in kpc
    Example: xGC, yGC, zGC = gal2XYZ(gl, gb, helio_dist)
    """
#    Rsun = 8.27 # in kpc (McMillan et al. 2011, Schoenrich et al. 2012)
#    Rsun = 8.00
    xGC = helio_dist*np.cos(np.radians(gb))*np.cos(np.radians(gl+180))+Rsun
    yGC = helio_dist*np.cos(np.radians(gb))*np.sin(np.radians(gl+180))
    zGC = helio_dist*np.sin(np.radians(gb))
    return xGC, yGC, zGC

def XYZ2gal(xGC, yGC, zGC, Rsun=8.00):
    """
    Description: Convert Galactocentric coordinates xyz into gl, gb, and heliocentric distance in kpc
    Usage: gl, gb, helio_dist = XYZ2gal(xGC, yGC, zGC)
    """
    ## distance to the Galactic center
    xnew = xGC - Rsun
    rcyl = np.sqrt(xnew**2+yGC**2)
    helio_dist = np.sqrt(xnew**2+yGC**2+zGC**2)
    gl = np.degrees(np.arctan2(yGC, xnew)) + 180.
    gl = np.where(gl >= 360, gl - 360., gl)
    if gl.size == 1:
        gl = np.atleast_1d(gl)[0]
    gb = np.degrees(np.arctan2(zGC, rcyl))
    return gl, gb, helio_dist

def eqPM2galPM(ra, dec, mu_ra, mu_dec):
    """
    Description: Given coordinates and proper motion in the equatorial J2000.0 system,
    return proper motions in the galactic coordinate system
    """

    # we need galactic coordinates
    gl, gb = eq2gal(ra, dec)
      
    ## angle phi between the NCP and the NGP (eqs.1-11)
    # but first, the equatorial coordinates of the NGP (J2000, from Wikipedia)
    alpha0 = 192.85948
    delta0 = 27.12825
    sinphi = np.sin(np.radians(90-dec))*np.sin(np.radians(ra-alpha0))/np.sin(np.radians(90-gb))
    aux1 = np.cos(np.radians(90-delta0))*np.sin(np.radians(90-dec))
    cosphi = (aux1-np.sin(np.radians(90-delta0))*np.cos(np.radians(90-dec))*np.cos(np.radians(ra-alpha0)))/np.sin(np.radians(90-gb))

    # rotation: eq.1-10
    mu_gl = mu_ra * cosphi + mu_dec * sinphi
    mu_gb = -mu_ra * sinphi + mu_dec * cosphi

    return mu_gl, mu_gb



def giFeH2Mr(gi, FeH):
        """
        Description: Calculate absolute magnitude in the SDSS r band for a given SDSS g-i color and FeH
        Example: Mr = giFeH2Mr(gi, FeH)
        """
        Mrgi = -0.56 + 14.32*gi - 12.97*gi**2. + 6.127*gi**3. - 1.267*gi**4. + 0.0967*gi**5. # TOMO II
        MrgiDelta = -1.11*FeH - 0.18*FeH**2. # TOMO II
        return Mrgi + MrgiDelta

def uggr2FeH(ug, gr):
        """
        Description: Calculate photometric metallicity from u-g and g-r colors using Bond et al.(2010) prescription
        Example: FeH = uggr2FeH(ug, gr)
        """
        return -13.13 + 14.09*ug + 28.04*gr -5.51*ug*gr -5.90*ug**2. -58.68*gr**2. + 9.14*(ug**2.)*gr -20.61*ug*gr**2. + 58.20*gr**3.

def gr2ri(gr):
        """
        Description: Interpolate r-i color from g-r color using Juric et al.(2008) prescription
        Example: ri = gr2ri(gr)
        """
        a = 4.9
        b = 2.45
        c = 1.68
        d = 0.050
        f = 1.39
        aa = (b**2.-3.*a*c)
        a0 = 2.*b**3. - 9.*a*b*c + 27.*a**2.*d + 27.*a**2.*np.log(1.-gr/f)
        a1 = np.sqrt(-4.*aa**3. + a0**2.)
        a2 = (-a0+a1)**(1./3.)
        a3 = -2.*b + 2.*2.**(1./3.) * aa / a2
        a4 = 2.**(2./3.)*a2
        return 1./(6.*a)*(a3 + a4)

def med_bin(x, y, min_x, max_x, binsize, return_empty=False):
    grid = np.arange(min_x, max_x, binsize)
    xC = grid + binsize/2.
    med = xC*0. - 99.99
    rms = xC*0. - 99.99
    cnt = xC*0. - 99.99
    medErr = xC*0. - 99.99
    for i in np.arange(xC.size):
        bin = np.where((x > grid[i]) & (x <= grid[i]+binsize))
        count = len(bin[0])
        if count > 4:
            med[i] = np.median(y[bin])
            rms[i] = 0.741*(np.percentile(y[bin],75) - np.percentile(y[bin],25))
            medErr[i] = np.sqrt(np.pi/2)*rms[i]/np.sqrt(count-1)
        cnt[i] = count
    if return_empty:
        return xC, med, rms, medErr, cnt
    else:
        index = np.where( cnt > -99 )
        xC = xC[index]
        med = med[index]
        rms = rms[index]
        medErr = medErr[index]
        cnt = cnt[index]
        return xC, med, rms, medErr, cnt

def systematic_velocity(rv, phase, ampR, phaseL, phaseR, fudge_factor=1.64):
    """ Description: Given measured heliocentric radial velocity,
        phase of pulsation and V band amplitude, provide the systemic
        radial velocity (center-of-mass velocity) using 3 different methods

        Usage: systematic_velocity(helio_rv, phase, ampV, phaseL, phaseR)
    """
    # get systemic velocity using X Ari template from Layden (1994)
    if phase < 0.77:
        k = 119.5
        l = rv - k*phase
        rvErr_layden = k*np.abs(phaseL - phaseR)
        rv_sys_layden = k*0.5+l
    else:
        k = -13.0
        rvErr_layden = k*np.abs(phaseL - phaseR)
        l = rv - k*phase
        y = k*0.77 + l
        k = 119.5
        l = y - k*0.77
        rv_sys_layden = k*0.5+l

    # get systemic velocity using X Ari template from Sesar et al.(2010b)
    k = 119.091
    l = rv - k*phase
    rv_sys_ses10b = k*0.5 + l
    rvErr_ses10b = k*np.abs(phaseL - phaseR)

    # get systemic velocity using slope from Kollmeier et al.(2009)
    k = 92.7
    l = rv - k*phase
    rv_sys_koll09 = k*0.5 + l
    rvErr_koll09 = k*np.abs(phaseL - phaseR)

    ## get systemic velocity using Liu (1991) synthetic light curve
    # load light curve data
    phase_template = np.arange(50)*0.02
    rv_template = np.genfromtxt('/home/bsesar/opt/python/synthetic_rv.dat', usecols=(1))
#    rv_template -= rv_template.min()
#    rv_template /= rv_template.max()
    phase_s = phase_template - 1
    phase_s = np.append(phase_s, phase_template)
    phase_s = np.append(phase_s, phase_template + 1)
#    ampV = 1.17*ampR # +- 0.05, from Sesar et al.(2012)
    ampV = 1.25*ampR # +- 0.05, from Sesar et al.(2012)
#    ampV = ampR
#    amp_RV = (40.5*ampV + 42.7)/0.89 # correlation between V band amplitude and RV amplitude from Liu (1991)
    amp_RV = 1.47*(40.5*ampV + 42.7)/1.35 # correlation between V band amplitude and RV amplitude from Liu (1991)
    rv_template *= amp_RV
    mag_s = np.append(np.append(rv_template,rv_template),rv_template)
    s = scipy.interpolate.UnivariateSpline(phase_s, mag_s, k=3, s=0)
    offset = rv - s(phase)
#    rv_sys_synth = s(0.553) + offset
    rv_sys_synth = s(0.5) + offset
#    rv_sys_synth = s(0.35) + offset
    rvErr_synth = np.abs(s(phaseL) - s(phaseR))
#    return rv_sys_layden[0], rv_sys_synth[0], rv_sys_koll09[0], rvErr_layden, rvErr_synth, rvErr_koll09
    return rv_sys_layden[0], rvErr_layden, rv_sys_synth[0], rvErr_synth

def layden_FeH(K, Hd, Hg, Hb, Keck=False):
    a = 13.858
    b = -1.185
    c = 4.228
    d = -0.32
    if Keck:
        K = 1.16*K - 0.39 # k = 1.16 +- 0.04, l = -0.39 +- 0.25, rms = 0.25
        Hd = 1.02*Hd + 1.39 # k = 1.02 +- 0.05, l = 1.39 +- 0.28, rms = 0.3
        Hg = 1.07*Hg - 0.45 # k = 1.07 +- 0.04, l = -0.45 +- 0.28, rms = 0.27
        Hb = 1.01*Hb + 1.18 # k = 1.01 +- 0.03, l = 1.18 +- 0.21, rms = 0.14
    else:
        K = 1.23*K - 0.44
        Hd = 0.84*Hd + 1.64
        Hg = 1.13*Hg - 0.63
        Hb = 1.02*Hb + 0.98
    WK = K
    WH3 = np.average([Hd, Hg, Hb])
    FeH_H3z = (WK - a - b*WH3)/(c+d*WH3)

    a = 13.669
    b = -1.125
    c = 4.147
    d = -0.3
    WH2 = np.average([Hd, Hg])
    FeH_H2 = (WK - a - b*WH2)/(c+d*WH2)

    return FeH_H3z, FeH_H2, Hd, Hg, Hb, WH3, WH2

def systematic_velocity2(rv_gamma, rvErr_gamma, rv_beta, rvErr_beta, phase, phaseL, phaseR, ampR, fudge_factor=1.21):
    """ Description: Given measured heliocentric radial velocity,
        phase of pulsation and V band amplitude, provide the systemic
        radial velocity (center-of-mass velocity) using 3 different methods

        Usage: systematic_velocity(helio_rv, phase, ampV, phaseL, phaseR)
    """
    ## get systemic velocities using Sesar (2012) synthetic RV curves
    ampV = fudge_factor*ampR
    amp_RV_gamma = 46.1*ampV + 38.5 # from Sesar (2012)
    amp_RV_beta = 42.1*ampV + 51.1 # from Sesar (2012)
    # load H-gamma RV template
    gamma_template = np.genfromtxt('/home/bsesar/projects/rrlyr/Hgamma.txt', names='phase, rv', dtype='f4,f4')
    gamma_template['rv'] *= amp_RV_gamma
    s_gamma = scipy.interpolate.UnivariateSpline(gamma_template['phase'], gamma_template['rv'], k=1, s=0)
#    offset_gamma = rv_gamma - s_gamma(phase)
    phase_grid = np.arange(phaseL, phaseR, 0.00001)
    if phaseL-phaseR > 0.5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    diff_gamma = np.abs(simps(s_gamma(phase_grid), phase_grid)/np.abs(phaseR-phaseL) - s_gamma(phase))
    offset_gamma = rv_gamma - (simps(s_gamma(phase_grid), phase_grid)/np.abs(phaseR-phaseL))
    rv_sys_gamma = gamma_template['rv'][135] + offset_gamma # phase_sys = 0.27
    #rv_sys_gamma = gamma_template['rv'][175] + offset_gamma # phase_sys = 0.35
    #rv_sys_gamma = gamma_template['rv'][250] + offset_gamma # phase_sys = 0.5
    rvErr_sys_gamma = (amp_RV_gamma+2.8)*np.sqrt(0.04**2 + (0.1*1.42)**2)
    rvErr_tot_gamma = np.sqrt(rvErr_gamma**2 + rvErr_sys_gamma**2)
    # load H-beta RV template
    beta_template = np.genfromtxt('/home/bsesar/projects/rrlyr/Hbeta.txt', names='phase, rv', dtype='f4,f4')
    beta_template['rv'] *= amp_RV_beta
    s_beta = scipy.interpolate.UnivariateSpline(beta_template['phase'], beta_template['rv'], k=1, s=0)
#    offset_beta = rv_beta - s_beta(phase)
    diff_beta = np.abs(simps(s_beta(phase_grid), phase_grid)/np.abs(phaseR-phaseL) - s_beta(phase))
    if (diff_gamma > 1) or (diff_beta > 1):
        print "Discrepant velocities!"
        sys.exit(1)
    offset_beta = rv_beta - (simps(s_beta(phase_grid), phase_grid)/np.abs(phaseR-phaseL)) 
    rv_sys_beta = beta_template['rv'][135] + offset_beta # phase_sys = 0.27
    #rv_sys_beta = beta_template['rv'][175] + offset_beta # phase_sys = 0.35
    #rv_sys_beta = beta_template['rv'][250] + offset_beta # phase_sys = 0.5
    rvErr_sys_beta = (amp_RV_beta+3.0)*np.sqrt(0.04**2 + (0.1*1.54)**2)
    rvErr_tot_beta = np.sqrt(rvErr_beta**2 + rvErr_sys_beta**2)
    rv_mean, sum_weights = np.average(np.array([rv_sys_gamma, rv_sys_beta]), weights=1./np.array([rvErr_tot_gamma, rvErr_tot_beta])**2, returned=True)
    rv_meanErr = 1./np.sqrt(sum_weights)
    spec_blurring_gamma = np.abs(s_gamma(phaseL) - s_gamma(phaseR))
    spec_blurring_beta = np.abs(s_beta(phaseL) - s_beta(phaseR))
    return rv_sys_gamma, rvErr_sys_gamma, rv_sys_beta, rvErr_sys_beta, abs(rv_sys_gamma-rv_sys_beta), rv_mean, rv_meanErr, spec_blurring_gamma, spec_blurring_beta, s_gamma(phase)-s_beta(phase)

def systematic_velocity_alpha(rv_alpha, rvErr_alpha, phase, phaseL, phaseR, ampR, scaling_factor=1.21):
    """ Description: Given measured heliocentric radial velocity,
        phase of pulsation and R band amplitude, provide the systemic
        radial velocity (center-of-mass velocity) using Sesar (2012) RV templates

        Usage: systematic_velocity_alpha(helio_rv, phase, ampR, phaseL, phaseR)
    """
    ## get systemic velocities using Sesar (2012) synthetic RV curves
    ampV = scaling_factor*ampR
    amp_RV_alpha = 35.6*ampV + 78.2 # from Sesar (2012)
    # load H-alpha RV template
    alpha_template = np.genfromtxt('/home/bsesar/papers/RV_templates/Halpha.txt', names='phase, rv', dtype='f4,f4')
    alpha_template['rv'] *= amp_RV_alpha
    s_alpha = scipy.interpolate.UnivariateSpline(alpha_template['phase'], alpha_template['rv'], k=1, s=0)
    offset_alpha = rv_alpha - s_alpha(phase)
    rv_sys_alpha = alpha_template['rv'][135] + offset_alpha # phase_sys = 0.27
    #rv_sys_alpha = alpha_template['rv'][175] + offset_alpha # phase_sys = 0.35
    #rv_sys_alpha = alpha_template['rv'][250] + offset_alpha # phase_sys = 0.5
    rvErr_sys_alpha = (amp_RV_alpha+3.4)*np.sqrt(0.08**2 + (0.1*1.47)**2)
    rvErr_tot_alpha = np.sqrt(rvErr_alpha**2 + rvErr_sys_alpha**2)
    spec_blurring_alpha = np.abs(s_alpha(phaseL) - s_alpha(phaseR))
    return rv_sys_alpha, rvErr_sys_alpha, rvErr_tot_alpha, spec_blurring_alpha

def gal2Orphan(gl, gb):
    # define Orphan stream orbital plane (from Newberg et al. 2010)
    phi = np.radians(128.79)
    theta = np.radians(54.39)
    psi = np.radians(90.7)
    M = np.array([[np.cos(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.sin(psi), np.cos(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi),  np.sin(psi)*np.sin(theta)],
                  [-np.sin(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.cos(psi), -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi),  np.cos(psi)*np.sin(theta)],
                  [np.sin(theta)*np.sin(phi), -np.sin(theta)*np.cos(phi), np.cos(theta)]             
                 ])
    B = np.array([[np.cos(np.radians(gb))*np.cos(np.radians(gl))],
                  [np.cos(np.radians(gb))*np.sin(np.radians(gl))],
                  [np.sin(np.radians(gb))]
                 ])
    A = np.dot(M, B)
    l_Orphan = np.degrees(np.arctan(A[1]/A[0]))
    b_Orphan = np.degrees(np.arcsin(A[2]))
    # correct the equator of the Orphan stream for l_Orphan < -15
    if l_Orphan[0] < -15:
        b_Orphan[0] = b_Orphan[0] + 0.00628*l_Orphan[0]**2 + 0.42*l_Orphan[0] + 5
    return l_Orphan[0], b_Orphan[0]

def Orphan2gal(l_Orphan, b_Orphan):
    # correct the equator of the Orphan stream for l_Orphan < -15
    if l_Orphan < -15:
        b_Orphan = b_Orphan - 0.00628*l_Orphan**2 - 0.42*l_Orphan - 5
    # define Orphan stream orbital plane (from Newberg et al. 2010)
    phi = np.radians(128.79)
    theta = np.radians(54.39)
    psi = np.radians(90.7)
    M = np.array([[np.cos(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.sin(psi), np.cos(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi),  np.sin(psi)*np.sin(theta)],
                  [-np.sin(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.cos(psi), -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi),  np.cos(psi)*np.sin(theta)],
                  [np.sin(theta)*np.sin(phi), -np.sin(theta)*np.cos(phi), np.cos(theta)]
                 ])
    B = np.array([[np.cos(np.radians(b_Orphan))*np.cos(np.radians(l_Orphan))],
                  [np.cos(np.radians(b_Orphan))*np.sin(np.radians(l_Orphan))],
                  [np.sin(np.radians(b_Orphan))]
                 ])
    A = np.dot(np.linalg.inv(M), B)
    gl = np.degrees(np.arctan(A[1]/A[0]))+180
    gb = np.degrees(np.arcsin(A[2]))
    return gl[0], gb[0]

def eq2Sgr(ra, dec):
    ra = np.radians(ra)
    dec = np.radians(dec)
    lambda_Belokurov2013 = np.degrees(np.arctan2(-0.93595354*np.cos(ra)*np.cos(dec) - 0.31910658*np.sin(ra)*np.cos(dec) + 0.14886895*np.sin(dec), 0.21215555*np.cos(ra)*np.cos(dec) - 0.84846291*np.sin(ra)*np.cos(dec)-0.48487186*np.sin(dec)))
    beta_Belokurov2013 = np.degrees(np.arcsin(0.28103559*np.cos(ra)*np.cos(dec) - 0.42223415*np.sin(ra)*np.cos(dec) + 0.86182209*np.sin(dec)))
    return lambda_Belokurov2013, beta_Belokurov2013

def Sgr2eq(lambda_Belokurov2013, beta_Belokurov2013):
    lambda_Belokurov2013 = np.radians(lambda_Belokurov2013)
    beta_Belokurov2013 = np.radians(beta_Belokurov2013)
    ra = np.degrees(np.arctan2(-0.84846291*np.cos(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) - 0.31910658*np.sin(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) - 0.42223415*np.sin(beta_Belokurov2013), 0.21215555*np.cos(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) - 0.93595354*np.sin(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) + 0.28103559*np.sin(beta_Belokurov2013)))
    dec = np.degrees(np.arcsin(-0.48487186*np.cos(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) + 0.14886895*np.sin(lambda_Belokurov2013)*np.cos(beta_Belokurov2013) + 0.86182209*np.sin(beta_Belokurov2013)))
    return ra, dec

def lb_stream(ra, dec, ra_pole, dec_pole, ra0):
    node = 90 + ra_pole
    incl = 90 - dec_pole
    dr = np.radians(dec)
    angdiff = np.radians(ra - node)
    x1 = np.cos(angdiff)*np.cos(dr)
    y1 = np.sin(angdiff)*np.cos(dr)
    z1 = np.sin(dr)
    x2 = x1
    y2 = y1*np.cos(np.radians(incl)) + z1*np.sin(np.radians(incl))
    z2 =-y1*np.sin(np.radians(incl)) + z1*np.cos(np.radians(incl))
    b = np.degrees(np.arcsin(z2))
    l = np.degrees(np.arctan2(y2,x2)) + ra0
    return l, b

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
    print B11
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


def UVW2Vpmeq(U,V,W,ra,dec,D):
    '''#UVW2Vpmeq converts radial 3D velocities to velocities and proper motions.
    #Parameters:
    #   U, V, W: 3D velocity components with respect to LSR
    #   ra,dec:equatorial coordinates
    #   D: distance to the Sun in kpc.
    #Return:
    #   Vr: radial velocity from observation
    #   pmra, pmdec: proper motion in equatorial coordinates
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

    B21= 0.493076*cra*cdec-0.460007*sra*cdec+0.738424*sdec
    B22=-0.493076*sra-0.460007*cra+0
    B23=-0.493076*cra*sdec+0.460007*sra*sdec+0.738424*cdec
    B31=-0.867666*cra*cdec-0.198076*sra*cdec+0.455984*sdec
    B32= 0.867666*sra-0.198076*cra+0
    B33= 0.867666*cra*sdec +0.198076*sra*sdec +0.455984*cdec

    UVW =  np.array([U-velsun[0], V-velsun[1], W-velsun[2]])
    Vr = np.zeros(len(ra))
    pmra = np.zeros(len(ra))
    pmdec = np.zeros(len(ra))
    for idx in range(0, len(ra)):
	B = np.array([[B11[idx], k*D[idx]*B12[idx], k*D[idx]*B13[idx]], \
	    [B21[idx], k*D[idx]*B22[idx], k*D[idx]*B23[idx]], \
	    [B31[idx], k*D[idx]*B32[idx], k*D[idx]*B33[idx]]])
        iV = np.dot(inv(B), UVW)
	Vr[idx] = iV[0]
	pmra[idx] = iV[1]
	pmdec[idx] = iV[2]
    return Vr, pmra, pmdec

def bin2D(data,c1,c2,rmin,rmax,cmin,cmax,rstep,cstep):
    #data col1:row, col2:column
    #rmin,rmax:row range of binning matrix
    #cmin,cmax:column range of binning matrix
    #rstep,cstep:row/column bin width
    #length:data size
    length=len(data)
    dr=rmax-rmin
    dc=cmax-cmin
    nr=round(dr/rstep)
    nc=round(dc/cstep)

    B=np.zeros((nc,nr))
    for i in range(0,length):
        ir=round((data[i][c1]-rmin)/rstep)+1
        ic=round((data[i][c2]-cmin)/cstep)+1
        if ((ir>=0) & (ir<nr) & (ic>=0) & (ic<nc)):
            B[ic][ir]=B[ic][ir]+1      
    return B
