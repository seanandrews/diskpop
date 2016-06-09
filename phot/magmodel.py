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

# parse dependent variable (passbands) (must be better way)
bands = ascii.read('bandinfo.dat')
xpb = (bands['name'])[np.in1d(bands['name'], x)]
xwl = (bands['wl'])[np.in1d(bands['name'], x)]

# parse parameters: 
#    theta = [teff, logg, rstar, zstar, pi, Av]
#                0,    1,     2,     3,  4,  5]

# calculate absolute magnitudes (Eq 1)
Mabs = lib_Mabs(xpb, p[0], p[1], p[3], p[2])
print(Mabs)

# calculate extinction values at each band
A_lambda = extinct(xwl, p[5])
print(A_lambda)

# now convert to observed apparent magnitudes (Eq 2)
Mapp = Mabs - 5.*np.log10(p[4]) - 5. + A_lambda
print(Mapp)

# now map these back onto the input bands (must be better way)
Mapp_back = np.zeros(len(x), dtype='float')
for ib in np.arange(len(Mapp)):
    Mapp_back[np.in1d(x, xpb[ib])] = Mapp[ib]

print(Mapp_back)
