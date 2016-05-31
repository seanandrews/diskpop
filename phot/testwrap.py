import numpy as np
import matplotlib.pyplot as plt
from load_phoenix_spec import load_phoenix_spec
from scipy.interpolate import interp1d

# load a spectrum
filename = 'lte580-4.0-0.0a+0.0.BT-Settl.7'
wl, flam = load_phoenix_spec(filename)

# filter profile
bwl, brlam = np.loadtxt('B.dat').T

# interpolate filter onto model spectrum
fint = interp1d(bwl, brlam, bounds_error=False, fill_value=0.)
br_int = fint(wl)

# calculate a broadband flux density
br_flux = np.trapz(flam * br_int * wl, wl) / np.trapz(br_int * 1e4

print(br_flux)


plt.axis([0.1, 1., 0., 1.])
plt.plot(bwl,brlam)
plt.plot(wl, br_int, '--r')
plt.show()
