#
#
#

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import sys


# effective temperature prior

# inputs
Sbar  = 60.
eSbar = 1.
Tinput = 8700.

# load spectral type |-> temperature conversion file
dt = {'ST': np.str, 'STix': np.float64, 'Teff':np.float64, 'eTeff':np.float64}
a  = pd.read_csv('data/adopted_spt-teff.txt', dtype=dt,
                 names=['ST','STix','Teff','eTeff'])

# discretized relationship
S_g  = np.array(a['STix'])
T_g  = np.array(a['Teff'])
eT_g = np.array(a['eTeff'])

# need to interpolate for appropriate integration
tint = interp1d(S_g, T_g)
eint = interp1d(S_g, eT_g)
S  = np.linspace(np.min(S_g), np.max(S_g), num=10.*len(S_g))
T  = tint(S)
eT = eint(S)

# calculate p(S)
p_S = np.exp(-0.5*((S-Sbar)/eSbar )**2) / (np.sqrt(2.*np.pi)*eSbar)

# now calculate p(T)
p_T = np.zeros_like(T)
for i in np.arange(len(T)): 
    p_TS = np.exp(-0.5*((T[i]-tint(S))/eint(S))**2) /  \
           (np.sqrt(2.*np.pi)*eint(S))
    p_T[i] = np.trapz(p_TS*p_S, S)

# create an interpolator for p_T
p_tint = interp1d(T, p_T)
prior_T = p_tint(Tinput)

print(prior_T)
