#-------------------------------------------------------------------------------
# The file defines the hydrostatic water pressure applied on the ice-ocean interface.
#-------------------------------------------------------------------------------

import numpy as np
from params import tides,tide_amplitude,tides_modulation

#---------------------Define sea level change timeseries------------------------
def sl_change(t):
    if tides == 'on':
    # Turn on the semi-diurnal tide (tidal frequency of 2 per day)
        if tides_modulation =='off':
            # Fixed tidal amplitude
            SLC = tide_amplitude*np.sin(4*np.pi*t/(3.154e7/12.0/30.0))  
        elif tides_modulation =='on':
            # Solar tidal modulation
            SLC = tide_amplitude/2.498*np.sin(2.0*np.pi*t/(0.5*3.154e7/12.0/30.0))\
            +1.5/2.498*tide_amplitude*np.sin(2.0*np.pi*t*12.42/12/(0.5*3.154e7/12.0/30.0)) 
    else:
        SLC = 0.0                                    # no sea level change for
                                                     # long-time marine problem
    return SLC
