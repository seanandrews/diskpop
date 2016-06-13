#
# Convert a magnitude to a flux density, or vice versa.
#

import numpy as np
import pandas as pd
import sys

def mag_to_flux(xb, xm, reverse=False, wl=False):

    # load information about available bands in library
    dt = {'name': np.str, 'wl': np.float64, 'zp': np.float64}
    a = pd.read_csv('bandinfo.txt', names=['name', 'wl', 'zp'], dtype=dt)
    band_nm = np.array(a['name'])
    band_wl = np.array(a['wl'])
    band_zp = np.array(a['zp'])

    # (slow) association of input bands with library bands and calculation of
    # flux densities or magnitudes
    if (reverse == False):
        flx = np.zeros_like(xm)
        owl = np.zeros_like(xm)
        for ix in np.arange(len(xb)):
            flx[ix] = band_zp[np.where(xb[ix] == band_nm)]*10.**(-0.4*xm[ix])
            owl[ix] = band_wl[np.where(xb[ix] == band_nm)]
        if (wl == False):
            return(flx)
        else: return( flx, owl )
    else:
        mag = np.zeros_like(xm)
        owl = np.zeros_like(xm)
        for ix in np.arange(len(xb)):
            mag[ix] = -2.5*np.log10(xm[ix]/band_zp[np.where(xb[ix] == band_nm)])
            owl[ix] = band_wl[np.where(xb[ix] == band_nm)]
        if (wl == False):
            return(mag)
        else: return( mag, owl )
