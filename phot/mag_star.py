# 
# Extract a list of model absolute magnitudes for a given {T, logg, Z, R} in 
# a specified set of bandpasses.  Then convert those into apparent magnitudes.
#

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from extinct import extinct

def mag_star(x, p):

    # parameter key: 
    #    p = [teff, logg, zstar, rstar, pi, Av]
    #            0,    1,     2,     3,  4,  5]

    # properly organize input bands (x), by sorting and uniqueing
    # save the shuffling, so you can re-format back to inputs in the end
    # this will also parse out an appropriate wavelength array for extinction
    #
    # load information about available bands in library
    dt = {'name': np.str, 'wl': np.float64, 'zp': np.float64}
    a = pd.read_csv('bandinfo.txt', names=['name', 'wl', 'zp'], dtype=dt)
    band_nm = np.array(a['name'])
    band_wl = np.array(a['wl'])
    #
    # sort the unique inputs and record that ordering for reshuffle later
    xs  = np.argsort(x)
    xus = np.argsort(xs)
    u_x, u_ix = np.unique(x[xs], return_inverse=True)
    bsel = np.in1d(np.array(band_nm), np.array(u_x))
    psel = np.argsort(band_nm[bsel])
    xin = band_nm[bsel]
    xwl = band_wl[bsel]
        
    # import the model library
    lib = np.load('maglib.npz')
    lteff  = lib['teff']
    llogg  = lib['logg']
    lzstar = lib['zstar']
    lmlib  = lib['mlib']
    lband  = lib['band']

    # extract the relevant bands
    maglib = lmlib[:,:,:,np.in1d(lband, xin)]

    # if the input parameters are exactly on the model library grid, then just
    # use the grid; if not, then do the trilinear interpolation
    if (np.any(lteff==p[0]) & np.any(llogg==p[1]) & np.any(lzstar==p[2])):
        mag_int = maglib[lteff==p[0], llogg==p[1], lzstar==p[2]]
    else:
        fint = RegularGridInterpolator((lteff, llogg, lzstar), maglib)
        mag_int = fint(np.array([p[0], p[1], p[2]]))

    # convert to standard absolute magnitudes
    Mabs = np.squeeze(mag_int) - 5.*np.log10(p[3]*6.96e10/3.0857e18) + 5.

    # calculate extinctions at each band
    A_lambda = extinct(xwl, p[5])

    # distance scaling for apparent magnitudes
    mapp = Mabs - 5.*np.log10(p[4]) - 5. + A_lambda

    # reshuffle and populate an output apparent mag array in same ordering as
    # requested input array
    fmapp = ((mapp[psel])[u_ix])[xus]

    # return the star contribution
    return(fmapp) 
