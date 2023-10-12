#-------------------------------------------------------------------------------
# The file defines the bed topohraphy and initial ice-sheet geometry
#-------------------------------------------------------------------------------

import numpy as np
from params import Lngth,Hght,bed_slope

# Generate Bed Topography
def bed(x):
    # Default marine ice sheet bed geometry
    Bed = 0.5 * Lngth * bed_slope - x * bed_slope
    return Bed

# Generate initial ice-water/ice-bed interface
def interface(x):
    Int = 0.5*(bed(x) + np.abs(bed(x)))
    return Int

# Generate initial ice-sheet top surface
def surface(x):
    Int = 0.5*(bed(x) + np.abs(bed(x))) + Hght
    return Int