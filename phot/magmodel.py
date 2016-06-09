#
# Generate a set of model magnitudes.
#

import numpy as np
from astropy.io import ascii
from lib_Mabs import lib_Mabs
from extinct import extinct

def magmodel(x, p):

    # parameter key: 
    #    p = [teff, logg, rstar, zstar, pi, Av]
    #            0,    1,     2,     3,  4,  5]

    # parse dependent variable (passbands) (must be better way)
    bands = ascii.read('bandinfo.dat')
    xpb = (bands['name'])[np.in1d(bands['name'], x)]
    xwl = (bands['wl'])[np.in1d(bands['name'], x)]

    # calculate absolute magnitudes (Eq 1)
    Mabs = lib_Mabs(xpb, p[0], p[1], p[3], p[2])

    # calculate extinction values at each band
    A_lambda = extinct(xwl, p[5])

    # now convert to observed apparent magnitudes (Eq 2)
    Mapp = Mabs - 5.*np.log10(p[4]) - 5. + A_lambda

    # now map these back onto the input bands (must be better way)
    Mapp_back = np.zeros(len(x), dtype='float')
    for ib in np.arange(len(Mapp)):
        Mapp_back[np.in1d(x, xpb[ib])] = Mapp[ib]

    return(Mapp_back)
