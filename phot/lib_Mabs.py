# 
# Extract a list of model absolute magnitudes for a given {T, logg, R, Z} in 
# a specified set of bandpasses.
#

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import sys

def lib_Mabs(band, teff, logg, zstar, rstar):

    # import the model library
    lib = np.load('maglib.npz')
    lteff  = lib['teff']
    llogg  = lib['logg']
    lzstar = lib['zstar']
    lmlib  = lib['mlib']
    lband  = lib['band']
 
    # extract the relevant bands
    maglib = lmlib[:,:,:,np.in1d(lband, band)]

    # if the input parameters are exactly on the model library grid, then just
    # use the grid; if not, then do the trilinear interpolation
    if (np.any(lteff==teff) & np.any(llogg==logg) & np.any(lzstar==zstar)):
        mag_int = maglib[lteff==teff, llogg==logg, lzstar==zstar]
    else:
        fint = RegularGridInterpolator((lteff, llogg, lzstar), maglib)
        mag_int = fint(np.array([teff, logg, zstar]))

    # put stellar radius in parsec units
    rstar *= 6.96e10/3.0857e18

    # convert to standard absolute magnitudes
    M = np.squeeze(mag_int) - 5.*np.log10(rstar) + 5.

    return(M)
