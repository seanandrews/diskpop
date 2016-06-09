#
# Generate a set of model magnitudes.
#

import numpy as np
from astropy.io import ascii
from lib_Mabs import lib_Mabs
from extinct import extinct

# temporary hard-coded function inputs
x = np.array(['Vj', 'Ic', 'Vj', 'K2m'])
p = np.array([4083., 4.1, 1.0, 0.0, 1./140., 0.2])

# parse dependent variable (passbands)
bands = ascii.read('bandinfo.dat')
xpb = (bands['name'])[np.in1d(bands['name'], x)]
xwl = (bands['wl'])[np.in1d(bands['name'], x)]

# parse parameters: 
#    theta = [teff, logg, rstar, zstar, pi, Av]
#                0,    1,     2,     3,  4,  5]

# calculate absolute magnitudes
Mabs = lib_Mabs(xpb, 4000., 4.0, 0.0, 1.0)
print(Mabs)
