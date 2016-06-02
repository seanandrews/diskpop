#
# Calculate A_lambda (extinction in magnitudes at a given wavelength).
#

import numpy as np
from scipy.interpolate import interp1d

def extinct(wl, Av):

    # load extinction curve(s) 
    ext  = np.loadtxt('ext_curves.dat', usecols=(0,1)).T
    ewl  = ext[0]
    AlAk = ext[1]

    # interpolating function
    fint = interp1d(ewl, AlAk)
    
    # calculate A_lambda, given A_V (denominator = A_V / A_K)
    A_lambda = Av * fint(wl) / fint(0.5488)

    return(A_lambda)
