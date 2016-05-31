import numpy as np
import matplotlib.pyplot as plt
import sys


# load the Johnson colors
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.JOHNSON'
mags_johnson = np.loadtxt(fname, comments='!', usecols=(0,1,2,3,4,5,6)).T
teff = mags_johnson[0,1:]
logg = mags_johnson[1,1:]
mh   = mags_johnson[2,1:]
ah   = mags_johnson[3,1:]
M_U  = mags_johnson[4,1:]
M_B  = mags_johnson[5,1:]
M_V  = mags_johnson[6,1:]

# take only solar metallicity models to start
noah = (ah == 0.0) & (mh == 0.0)
teff = teff[noah]
logg = logg[noah]
M_U  = M_U[noah]
M_B  = M_B[noah]
M_V  = M_V[noah]

# unique values of each parameter
uteff = np.unique(teff)
ulogg = np.unique(logg)

# master arrays
magU = np.full((len(uteff), len(ulogg)), np.nan)
magB = np.full((len(uteff), len(ulogg)), np.nan)
magV = np.full((len(uteff), len(ulogg)), np.nan)

# loop through the logg's and populate the magnitude arrays
for ig in np.arange(len(ulogg)):

    # at a fixed logg, note where you have magnitudes
    ind_step = (logg == ulogg[ig])
    fteff = teff[ind_step]
    loc = np.in1d(uteff, fteff)

    # populate master vectors
    magU[loc,ig] = M_U[ind_step]
    magB[loc,ig] = M_B[ind_step]
    magV[loc,ig] = M_V[ind_step]



