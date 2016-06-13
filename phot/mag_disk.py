#
# Calculate the disk contributions to the photometry, from both veiling and
# thermal dust emission, based on my prescription.
#

import numpy as np
from astropy.io import ascii

def mag_disk(x, p, wl=0):

    # parameter key: 
    #    p = [r_V^acc, psi, r_K^dust, T_s]
    #               0,   1,        2,   3]

    # properly organize input bands (x), by sorting and uniqueing
    # save the shuffling, so you can re-format back to inputs in the end
    # this will also parse out an appropriate wavelength array for extinction
    #
    # load information about available bands in library
    dt = {'name': np.str, 'wl': np.float64, 'zp': np.float64}
    a = pd.read_csv('bandinfo.txt', names=['name', 'wl', 'zp'], dtype=dt)
    band_nm = np.array(a['name'])
    band_wl = np.array(a['wl'])
    #
    # sort the unique inputs and record that ordering for reshuffle later
    xs  = np.argsort(x)
    xus = np.argsort(xs)
    u_x, u_ix = np.unique(x[xs], return_inverse=True)
    bsel = np.in1d(np.array(band_nm), np.array(u_x))
    psel = np.argsort(band_nm[bsel])
    xwl = band_wl[bsel]
    
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

    # reshuffle and populate an output apparent mag array in same ordering as
    # requested input array
    fr_lambda = ((r_lambda[psel])[u_ix])[xus]

    # return the disk contribution
    return(-2.5*np.log10(1.+fr_lambda))
