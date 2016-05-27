import numpy as np

def load_phoenix_spec(filename):

    # Load the spectrum (note here I've manually changed the file to 'E' format 
    # instead of 'D' format: that sucks but I'll figure out a fix later.  This
    # means of loading the spectrum is 10x faster than astropy's ascii.read.
    wlang, logflux = np.loadtxt(filename, usecols=(0,1)).T

    # Convert wavelengths from vacuum to air
    wlang = wlang / (1.+2.73512e-4+(131.4182/wlang**2)+(2.87249e8/wlang**4))

    # Inexplicably, the data are not sorted by wavelength.  Sort and convert to 
    # sensical units (microns and erg/s/cm**2/micron).
    s_wl = wlang.argsort()
    wlum = 1e-4 * wlang[s_wl]
    Flam = 1e4 * 10.**(logflux[s_wl]-8.)

    spec = wlum, Flam

    return spec
