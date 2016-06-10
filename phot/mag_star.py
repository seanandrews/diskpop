# 
# Extract a list of model absolute magnitudes for a given {T, logg, Z, R} in 
# a specified set of bandpasses.  Then convert those into apparent magnitudes.
#

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from astropy.io import ascii
from extinct import extinct

def mag_star(x, p):

    # parameter key: 
    #    p = [teff, logg, zstar, rstar, pi, Av, ..., ...]
    #            0,    1,     2,     3,  4,  5, ..., ...]

    # import the model library
    lib = np.load('maglib.npz')
    lteff  = lib['teff']
    llogg  = lib['logg']
    lzstar = lib['zstar']
    lmlib  = lib['mlib']
    lband  = lib['band']

    # extract the relevant bands
    maglib = lmlib[:,:,:,np.in1d(lband, x)]

    # if the input parameters are exactly on the model library grid, then just
    # use the grid; if not, then do the trilinear interpolation
    if (np.any(lteff==p[0]) & np.any(llogg==p[1]) & np.any(lzstar==p[2])):
        mag_int = maglib[lteff==p[0], llogg==p[1], lzstar==p[2]]
    else:
        fint = RegularGridInterpolator((lteff, llogg, lzstar), maglib)
        mag_int = fint(np.array([p[0], p[1], p[2]]))

    # convert to standard absolute magnitudes
    Mabs = np.squeeze(mag_int) - 5.*np.log10(p[3]*6.96e10/3.0857e18) + 5.

    # get wavelengths for each passband (must be a better way...)
    bands = ascii.read('bandinfo.dat')
    xwl = (bands['wl'])[np.in1d(bands['name'], x)]

    # calculate extinctions at each band
    A_lambda = extinct(xwl, p[5])

    # distance scaling for apparent magnitudes
    mapp = Mabs - 5.*np.log10(p[4]) - 5. + A_lambda

    return(mapp) 
