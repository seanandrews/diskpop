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
noah = (ah == 0.0) 
mh   = mh[noah]
teff = teff[noah]
logg = logg[noah]
M_U  = M_U[noah]
M_B  = M_B[noah]
M_V  = M_V[noah]

# unique values of each parameter
umh   = np.unique(mh)
uteff = np.unique(teff)
ulogg = np.unique(logg)

# master arrays
magU = np.full((len(uteff), len(ulogg), len(umh)), np.nan)
magB = np.full((len(uteff), len(ulogg), len(umh)), np.nan)
magV = np.full((len(uteff), len(ulogg), len(umh)), np.nan)

# loop through the mh's and logg's and populate the magnitude arrays
for im in np.arange(len(umh)):

    for ig in np.arange(len(ulogg)):

        # at a fixed logg and mh, note which teffs have magnitudes
        ind_step = (logg == ulogg[ig]) & (mh == umh[im])
        fteff = teff[ind_step]
        loc = np.in1d(uteff, fteff)

        # populate master vectors
        magU[loc, ig, im] = M_U[ind_step]
        magB[loc, ig, im] = M_B[ind_step]
        magV[loc, ig, im] = M_V[ind_step]


plt.plot(np.log10(teff[(mh == umh[3]) & (logg == ulogg[4])]), \
         M_U[(mh == umh[3]) & (logg == ulogg[4])], 'or')
plt.plot(np.log10(uteff), magU[:, 4, 3], '.b')
plt.show()


