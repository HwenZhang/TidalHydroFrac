#-------------------------------------------------------------------------------
# Define bed topography and initial ice-water interface functions.
# Note: Bed and ice-water interface should be equal on margins of the domain
# for the lake problem! They should be equal on one margin of the domain for the
# marine ice sheet problem (i.e., the grounded portion).
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from params import Lngth,Hght,model,model_setup,bed_slope

#-------------------- Generate Bed Topography-----------------------------------
def bed(x):
    #-------------------Default marine ice sheet bed geometry-------------------
    if model == 'marine':
        # Generate linear bed topography
        # 1. Lngth=4e4
        # slope = 2.5e-4
        # Bed = 0.25 * Lngth * slope - x * slope
        # 2. Lngth=2e4
        Bed = 0.5 * Lngth * bed_slope - x * bed_slope
    return Bed
#------------------Generate initial ice-water/ice-bed interface-----------------
def interface(x):
    if model == 'lake':
        Int = 0.5*(bed(x) - 5 + np.abs(bed(x) - (-5)))

        if model_setup == 'wedge_test':
            Int = 0.5*(bed(x) - 7.5 + np.abs(bed(x) - (-7.5)))

    elif model == 'marine':
        Int = 0.5*(bed(x) + np.abs(bed(x)))
    return Int
#-------------------------------------------------------------------------------
