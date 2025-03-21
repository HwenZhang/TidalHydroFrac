import numpy as np

def to_scientific_notation(number):
    scientific_str = "{:.0e}".format(number)
    scientific_str = scientific_str.replace('.', '_')
    return scientific_str

#-------------------------------------------------------------------------------
#----------------------------- MODEL OPTIONS -----------------------------------
#-------------------------------------------------------------------------------
# plotting flag
model_setup = 'tides_paper'

model = 'marine'
resultsname = 'results' # directory name

# For the marine ice sheet setup, turn the tidal cycle 'on' or 'off'.
tides = 'on'
tides_modulation = 'off'     # solar tidal modulation

# Turn 'on' or 'off' Newton convergence information:
print_convergence = 'on'

# Newtonian fluid when True, and Glen's Flow law when False
newtonian = True

# Save the intermediate time series
checkpoint = True

# Mesh resolution at the lower boundary
DX_s = 25.0                   # Mesh Size (in meters)
DX_h = 250.0                  # Mesh Size at the upper surface (in meters)

#-------------------------------------------------------------------------------
#-----------------------------MODEL PARAMETERS----------------------------------
#-------------------------------------------------------------------------------
# physical units:
# time - seconds
# space - meters
# pressure - pascals
# mass - kg

# Material parameters
A0 = 1.2e-26                       # Glen's law coefficient (ice softness, Pa^{-n}/s)
n = 3.0                            # Glen's law exponent
rm2 = 1 + 1.0/n - 2.0              # Exponent in variational forms: r-2
rm2_s = -2.0/3.0                   # Exponent in sliding law
B0 = A0**(-1/n)                    # Ice hardness (Pa s^{1/n})
B = (2**((n-1.0)/(2*n)))*B0        # "2*Viscosity" constant in weak form (Pa s^{1/n})
rho_i = 917.0                      # Density of ice (kg/m^3)
rho_w = 1027.0                     # Density of water (kg/m^3)
g = 9.81                           # Gravitational acceleration (m/s^2)
C = 1.2e7                          # Sliding law friction coefficient (Pa s^{1/n}/m)
mu = 0.30e9                        # Shear modulus (Pa)
eta_const = 1.277e+15              # Newtonian viscosity if applied (Pa/s)

# Numerical parameters
eps_v = 1.0e-18                    # Flow law regularization parameter
eps_p = 1.0e-14                    # Penalty method parameter
eps_vs= 1.0e-15                    # Sliding law regularisation parameter
quad_degree = 16                   # Quadrature degree for weak forms

tol = 1.0e-3                       # Numerical tolerance for boundary geometry:
                                   # s(x,t) - b(x) > tol on ice-water boundary,
                                   # s(x,t) - b(x) <= tol on ice-bed boundary.
# Geometry parameters
Lngth = 20*1000.0                  # Length of the domain (m)
Hght = 500.0                       # (Initial) Height of the domain (m)
bed_slope = 2e-2                 # bedrock slope

tide_amplitude = 1.00              # tide amplitude (m)
sea_level = Hght*(917.0/1027.0)    # Sea level elevation (m).
                                   # (Initial sea level for the tides problem)     
                                                              
# Time-stepping parameters for tidal problems
if tides == 'off': # mesh initiation
    nt_per_year = 0.05 * 1e3           # mesh initiation time step
    t_final = 50 * 3.154e7             # 40 years
else: # tidal simulation
    if tides_modulation == 'on':
        nt_per_year =  50*1000             # Number of timesteps per year. 
        t_final =  40.0/360.*3.154e7       # Final time (yr*sec_per_year).
    else:
        nt_per_year =  50*1000             # Number of timesteps per year. 
        t_final =  20.0/360.*3.154e7       # Final time (yr*sec_per_year).

nt = int(nt_per_year*t_final/3.154e7) # Number of time steps
dt = t_final/nt                       # Timestep size

nx = 10*int(Lngth/DX_s)               # Horizontal coordinate for computing surface
X_fine = np.linspace(0,Lngth,num=nx)  # slopes, interpolated grounding line positions, and plotting.
save_interval = 100

# Set positive inflow speed boundary condition for marine ice sheet problem
U0 = 17.45/3.154e7                    # Inflow speed 1 km/yr (m/yr / sec/yr) ~ 3e-3

slope_str = 'Slope2e_2'
#------------------------------Casename----------------------------------------

if tides == 'off':
    if newtonian:
        casename = 'stokes_mesh_initiation_U'+"%02d" % int(U0*3.154e7)+'ma_L20000_'+slope_str+'_eta1_2e12_n3_0_'+\
                'mu0_30e9_deltap1e_14_deltav1e_18_C1_2e7_DX'+str(int(DX_s))
    else:
        casename = 'stokes_mesh_initiation_U'+"%02d" % int(U0*3.154e7)+'ma_L20000_'+slope_str+'_A1e_25_n3_0_'+\
                'mu0_30e9_deltap1e_13_deltav1e_18_C1_2e7_DX'+str(int(DX_s))
else:
    if tides_modulation == 'on':
        casename = 'stokes_tidal_response_U' + "%02d" % int(U0*3.154e7)+'ma_L20000_'+slope_str+'_A1e_25_n3_0_'+\
                    'mu0_30e9_deltap1e_13_deltav1e_18_C1_2e7_DX'+str(int(DX_s))+'_modtide' + ("%0.2f" % tide_amplitude).replace('.','_')
    elif newtonian:
        casename = 'stokes_tidal_response_U' + "%02d" % int(U0*3.154e7)+'ma_L20000_'+slope_str+'_eta1_2e15_n3_0_'+\
                    'mu0_30e9_deltap1e_14_deltav1e_18_C1_2e7_DX'+str(int(DX_s))+'_tide' + ("%0.2f" % tide_amplitude).replace('.','_')
    else:
        casename = 'stokes_tidal_response_U' + "%02d" % int(U0*3.154e7)+'ma_L20000_'+slope_str+'_A1e_25_n3_0_'+\
                    'mu0_30e9_deltap1e_14_deltav1e_18_C1_2e7_DX'+str(int(DX_s))+'_tide' + ("%0.2f" % tide_amplitude).replace('.','_')