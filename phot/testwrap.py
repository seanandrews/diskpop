import numpy as np
from load_phoenix_spec import load_phoenix_spec

filename = 'lte580-4.0-0.0a+0.0.BT-Settl.7'
wl, flam = load_phoenix_spec(filename)

print(wl)


