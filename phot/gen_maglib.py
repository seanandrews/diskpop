import numpy as np


# load the Johnson colors
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.JOHNSON.Vega'
mags_johnson = np.loadtxt(fname, comments='!', usecols=(0,1,2,3,4,5,6)).T
teff = mags_johnson[0,1:]
logg = mags_johnson[1,1:]
mh   = mags_johnson[2,1:]
ah   = mags_johnson[3,1:]
M_U  = mags_johnson[4,1:]
M_B  = mags_johnson[5,1:]
M_V  = mags_johnson[6,1:]

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

# Extract analogous Cousins RI data
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.COUSINS.Vega'
mags_cousins= np.loadtxt(fname, comments='!', usecols=(6,7)).T
M_R = (mags_cousins[0,1:])[noah]
M_I = (mags_cousins[1,1:])[noah]

# Extract analogous 2MASS JHKs data
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.2MASS.Vega'
mags_twomass= np.loadtxt(fname, comments='!', usecols=(4,5,6)).T
M_J = (mags_twomass[0,1:])[noah]
M_H = (mags_twomass[1,1:])[noah]
M_K = (mags_twomass[2,1:])[noah]

# Extract analogous Spitzer IRAC, MIPS-24 data
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.SPITZER.Vega'
mags_spitzer = np.loadtxt(fname, comments='!', usecols=(4,5,6,7,10)).T
M_IRAC1 = (mags_spitzer[0,1:])[noah]
M_IRAC2 = (mags_spitzer[1,1:])[noah]
M_IRAC3 = (mags_spitzer[2,1:])[noah]
M_IRAC4 = (mags_spitzer[3,1:])[noah]
M_MIPS1 = (mags_spitzer[4,1:])[noah]

# Extract analogous WISE data
fname = 'BT-Settl/AGS2009/colmag.BT-Settl.server.WISE.Vega'
mags_wise= np.loadtxt(fname, comments='!', usecols=(4,5,6,7)).T
M_WISE1 = (mags_wise[0,1:])[noah]
M_WISE2 = (mags_wise[1,1:])[noah]
M_WISE3 = (mags_wise[2,1:])[noah]
M_WISE4 = (mags_wise[3,1:])[noah]

# master arrays
magU = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magB = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magV = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magR = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magI = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magJ = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magH = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magK = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magI1 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magI2 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magI3 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magI4 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magM1 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magW1 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magW2 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magW3 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)
magW4 = np.full((len(uteff), len(ulogg), len(umh)), np.inf)


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
        magR[loc, ig, im] = M_R[ind_step]
        magI[loc, ig, im] = M_I[ind_step]
        magJ[loc, ig, im] = M_J[ind_step]
        magH[loc, ig, im] = M_H[ind_step]
        magK[loc, ig, im] = M_K[ind_step]
        magI1[loc, ig, im] = M_IRAC1[ind_step]
        magI2[loc, ig, im] = M_IRAC2[ind_step]
        magI3[loc, ig, im] = M_IRAC3[ind_step]
        magI4[loc, ig, im] = M_IRAC4[ind_step]
        magM1[loc, ig, im] = M_MIPS1[ind_step]
        magW1[loc, ig, im] = M_WISE1[ind_step]
        magW2[loc, ig, im] = M_WISE2[ind_step]
        magW3[loc, ig, im] = M_WISE3[ind_step]
        magW4[loc, ig, im] = M_WISE4[ind_step]


# now stack them!
magbar = np.stack((magU, magB, magV, magR, magI, magJ, magH, magK, magI1, magI2, magI3, magI4, magM1, magW1, magW2, magW3, magW4), axis=-1)


# package the results
filts = np.array(['Uj','Bj','Vj','Rc','Ic','J2m','H2m','K2m','IRAC1','IRAC2','IRAC3','IRAC4','MIPS1','WISE1','WISE2','WISE3','WISE4'])

wleff = np.array([0.3597, 0.4377, 0.5488, 0.6515, 0.7981, 1.235, 1.662, 2.159, 3.5508, 4.4960, 5.7245, 7.8843, 23.67, 3.3526, 4.6028, 11.5608, 22.0883])

fnu_zp = np.array([1801.4, 4092.6, 3698.3, 3124.6, 2502.7, 1594., 1024., 666.8, 280.9, 179.7, 115.0, 64.13, 7.14, 309.450, 171.787, 31.674, 8.363])

np.savez('maglib', teff=uteff, logg=ulogg, zstar=umh, band=filts, band_wl=wleff, fnu_zp=fnu_zp, mlib=magbar)
