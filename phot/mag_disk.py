#
# Calculate the disk contributions to the photometry, from both veiling and
# thermal dust emission, based on my prescription.
#

import numpy as np
from astropy.io import ascii

def mag_disk(x, p):

    # parameter key: 
    #    p = [r_V^acc, psi, r_K^dust, T_s]
    #               0,   1,        2,   3]

    # get wavelengths for each passband (must be a better way...)
    bands = ascii.read('bandinfo.dat')
    xwl = (bands['wl'])[np.in1d(bands['name'], x)]

    # calculate r_lambda^acc (Eq 4)
    Vwl = 0.55
    Jwl = 1.235
    bet = 10.
    r_acc = p[0] * (xwl/Vwl)**(-p[1]) * (1.- (1./(1.+np.exp(-bet*(xwl-Jwl)))))    
    # calculate r_lambda^dust (Eq 5)
    h = 6.626e-27
    c = 2.9979e14	# in microns/s
    k = 1.381e-16
    Kwl = 2.159
    r_dust = p[2] * (Kwl/xwl)**5 * (np.exp(h*c/(Kwl*k*p[3]))-1.) / \
                                   (np.exp(h*c/(xwl*k*p[3]))-1.)

    # combine the contributions
    r_lambda = r_acc + r_dust

    # return the disk contribution
    return(-2.5*np.log10(1.+r_lambda))
