# 
# Extract a list of model absolute magnitudes for a given {T, logg, R, Z} in 
# a specified set of bandpasses.
#


import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import sys

# parse inputs
in_band  = np.array(['Rc', 'Ic', 'J2m'])
in_teff  = 4800.
in_logg  = 4.5
in_zstar = 0.0
in_rstar = 1.0

# import the model library
lib = np.load('maglib.npz')
teff  = lib['teff']
logg  = lib['logg']
zstar = lib['zstar']
mlib  = lib['mlib']
band  = lib['band']

# extract the relevant bands
maglib = mlib[:,:,:,np.in1d(band, in_band)]

# if the input parameters are *exactly* on the model library grid, then just
# use the grid; if not, then do the interpolation
if (np.any(teff == in_teff) & \
    np.any(logg == in_logg) & \
    np.any(zstar == in_zstar)):

    mag_int = np.squeeze(maglib[teff==in_teff, logg==in_logg, zstar==in_zstar])

else:

    # create the interpolator function (here using trilinear interpolation)
    fint = RegularGridInterpolator((teff, logg, zstar), maglib)

    # calculate interpolated absolute magnitudes at stellar surface
    mag_int = np.squeeze(fint(np.array([in_teff, in_logg, in_zstar])))

print(mag_int)

# put stellar radius in parsec units
in_rstar *= 6.96e10/3.0857e18

# convert to standard absolute magnitudes
M = mag_int - 5.*np.log10(in_rstar) + 5.

print(M)
