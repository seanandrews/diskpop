# 
# Description...
#

import numpy as np
import pandas as pd
import time
from scipy.interpolate import RegularGridInterpolator as RGI
from extinct import extinct
import matplotlib.pyplot as plt
import sys

def mag_model(x, p, lib=None):

    # parameter key: 
    #    p = [teff, logg, zstar, rstar, pi, Av, r_V^acc, psi, r_K^dust, T_s]
    #            0,    1,     2,     3,  4,  5,       6,   7,        8,   9]

    # properly organize input bands (x), by sorting and uniqueing
    # save the shuffling, so you can re-format back to inputs in the end
    # this will also parse out an appropriate wavelength array for extinction
    #
    # load information about available bands in library
    dt = {'name': np.str, 'wl': np.float64, 'zp': np.float64}
    a = pd.read_csv('bandinfo.txt', names=['name', 'wl', 'zp'], dtype=dt)
    band_nm = np.array(a['name'])
    band_zp = np.array(a['zp'])
    band_wl = np.array(a['wl'])
    #
    # sort the unique inputs and record that ordering for reshuffle later
    xs  = np.argsort(x)
    xus = np.argsort(xs)
    u_x, u_ix = np.unique(x[xs], return_inverse=True)
    bsel = np.in1d(np.array(band_nm), np.array(u_x))
    psel = np.argsort(band_nm[bsel])
    xin = band_nm[bsel]
    xzp = band_zp[bsel]
    xwl = band_wl[bsel]
        
    # import the model library if you didn't pass it already
    if (lib is None): 
        lib = np.load('maglib.npz')		# bottleneck
        lteff  = lib['teff']
        llogg  = lib['logg']
        lzstar = lib['zstar']
        lmlib  = lib['mlib']
        lband  = lib['band']
    else: lteff, llogg, lzstar, lmlib, lband = lib

    # extract the relevant bands
    maglib = lmlib[:,:,:,np.in1d(lband, xin)]

    # if the input parameters are exactly on the model library grid, then just
    # use the grid; if not, then do the trilinear interpolation
    if (np.any(lteff==p[0]) & np.any(llogg==p[1]) & np.any(lzstar==p[2])):
        mag_int = maglib[lteff==p[0], llogg==p[1], lzstar==p[2]]
    else:
        fint = RGI((lteff, llogg, lzstar), maglib)
        mag_int = fint(np.array([p[0], p[1], p[2]]))

    # convert to standard absolute magnitudes
    Mabs = np.squeeze(mag_int) - 5.*np.log10(p[3]*6.96e10/3.0857e18) + 5.

    # apparent magnitudes --> flux densities (no extinction yet)
    fstar = xzp*10.**(-0.4*(Mabs - 5.*np.log10(p[4]) - 5.)) 

    # contribution from veiling
    wlV = np.squeeze(band_wl[np.where(band_nm=='Vj')])
    zpV = np.squeeze(band_zp[np.where(band_nm=='Vj')])
    wlJ = np.squeeze(band_wl[np.where(band_nm=='J2m')])
    bet = 10.
    if (np.any(xin == 'Vj')):
        fstarV = np.squeeze(fstar[np.where(xin=='Vj')])
        MabsV = np.squeeze(Mabs[np.where(xin=='Vj')])
    else:
        fVint = RGI((lteff, llogg, lzstar), lmlib[:,:,:,np.where(lband=='Vj')])
        MabsV = np.squeeze(fVint(np.array([p[0], p[1], p[2]])) - \
                5.*np.log10(p[3]*6.96e10/3.0857e18) + 5.)
        fstarV = np.squeeze(zpV*10.**(-0.4*(MabsV - 5.*np.log10(p[4])-5.)))
    #fveil = p[6]*fstarV + p[7]*(xwl-wlV)	# linear model
    fveil = p[6]*fstarV*(xwl/wlV)**p[7]		# power-law model
    fveil *= (1. - 1./(1.+np.exp(-bet*(xwl-wlJ))))	# taper veiling
    
    # convert to composite apparent magnitudes and redden appropriately
    A_lambda = extinct(xwl, p[5])
    mtot = -2.5*np.log10((fstar + fveil) / xzp) + A_lambda

    # reshuffle and populate an output apparent mag array in same ordering as
    # requested input array
    fmtot = ((mtot[psel])[u_ix])[xus]

    # return the star contribution
    return(fmtot) 
