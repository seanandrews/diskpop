# 
# Extract a list of model absolute magnitudes for a given {T, logg, R, Z} in 
# a specified set of bandpasses.
#


import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import sys

# parse inputs
in_band  = 'Rc'
in_teff  = 4800.
in_logg  = 4.5
in_zstar = 0.0

# import the model library
lib = np.load('maglib.npz')
teff  = lib['teff']
logg  = lib['logg']
zstar = lib['zstar']
mlib  = lib['mlib']
band  = lib['band']


# temp: pluck out the in_band
maglib = np.squeeze(mlib[:,:,:,(band == in_band)])


# interpolator
fint = RegularGridInterpolator((teff, logg, zstar), maglib)


pts = np.array([in_teff, in_logg, in_zstar])
MM = fint(pts)
print(np.inf-4.)
print(logg)
print(zstar)

#plt.axis([4000, 5000, -40, -35])
plt.axis([-0.5, 6., -40, -35])
#plt.plot(teff, maglib[:,(logg == 4.0),(zstar == -0.5)], 'ob')
print(np.shape(maglib[(teff == 4800.),:,(zstar == 0.0)]))
plt.plot(logg, np.squeeze(maglib[(teff == 4800.),:,(zstar == 0.0)]), 'og')
#plt.plot(in_logg, MM, 'or')

plt.show()
