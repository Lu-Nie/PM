import numpy as np

def ds2tp(ra, dec, raz, decz):
    """Projection of spherical coordinates onto tangent plane:
       "gnomonic" projection - "standard coordinates" (double precision)

       Given:
       RA,DEC spherical coordinates of point to be projected (in degrees)
       RAZ,DECZ spherical coordinates of tangent point (in degrees)

       Returned:
       XI,ETA rectangular coordinates on tangent plane (in degrees)
       J int status: 0 = OK, star on tangent plane
                     1 = error
    """
    sdecz = np.sin(decz*np.pi/180)
    sdec = np.sin(dec*np.pi/180)
    cdecz = np.cos(decz*np.pi/180)
    cdec = np.cos(dec*np.pi/180)
    radif = (ra-raz)*np.pi/180
    sradif = np.sin(radif)
    cradif = np.cos(radif)

    denom = sdec*sdecz + cdec*cdecz*cradif
    j = np.where(denom > 1e-06, 0, 1)

    xi = (cdec*sradif/denom)*180/np.pi
    eta = ((sdec*cdecz - cdec*sdecz*cradif)/denom)*180/np.pi
    #print "min(ra), max(ra), min(dec), max(dec):", min(ra), max(ra), min(dec), max(dec)
    #print "raz, decz:", raz, decz
    #print "min(xi), max(xi), min(eta), max(eta):", min(xi), max(xi), min(eta), max(eta)

    return xi, eta, j

def dranrm(angle):
    aux = np.mod(angle, 2*np.pi)
    aux = np.where(aux < 0, aux + 2*np.pi, aux)
    return aux

def dtp2s(xi, eta, raz, decz):
    """Transform tangent plane coordinates into spherical
       (double precision)

       Given:
       XI,ETA tangent plane rectangular coordinates (in radians)
       RAZ,DECZ spherical coordinates of tangent point (in degrees)

       Returned:
       RA,DEC spherical coordinates (in degrees)

       Called: dranrm
    """
    sdecz = np.sin(decz*np.pi/180.)
    cdecz = np.cos(decz*np.pi/180.)

    denom = cdecz - eta*sdecz

    ra = dranrm(np.arctan2(xi, denom) + raz*np.pi/180.)*180/np.pi
    dec = np.arctan2(sdecz + eta*cdecz, np.sqrt(xi**2 + denom**2))*180/np.pi

    return ra, dec

