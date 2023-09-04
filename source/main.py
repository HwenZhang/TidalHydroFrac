#-------------------------------------------------------------------------------
# This program solves a nonlinear Stokes flow problem describing groundning line
# dynamics in two settings:
#
# (1) A marine ice sheet, OR ...
# (2) A subglacial lake that is (potentially) undergoing water volume changes.
#
# Choose model setup (1) or (2) in the params.py fßile.
#
# See the JFM manuscript for a complete formulation of the problem.
#-------------------------------------------------------------------------------
import sys
sys.path.insert(0, './scripts')

from dolfin import *
from fenics import *
import matplotlib.pyplot as plt
import numpy as np
from stokes import stokes_solve_marine,get_zero_m
from geometry import interface,bed,surface
from meshfcns import mesh_routine
from plotting import *
import scipy.integrate as scpint
import os
from params import (rho_i,g,tol,t_final,nt_per_year,Lngth,Hght,nt,dt,model,
                    print_convergence,X_fine,nx,tides,DX_s,model_setup,test)ß
from params import (rho_i,g,tol,B,rm2,rho_w,C,eps_p,eps_v,sea_level,dt,
                    quad_degree,Lngth,U0,mu,tide_amplitude,resultsname,casename)

#--------------------Initial conditions-----------------------------------------
# Compute initial mean elevation of ice-water interface and initial lake volume.
s_mean0 = np.mean(interface(X_fine)[interface(X_fine)-bed(X_fine)>tol])
lake_vol_0 = scpint.quad(lambda x: interface(x)-bed(x),0,Lngth,full_output=1)[0]
#-------------------------------------------------------------------------------

path = []
path.append(resultsname+'/' + casename )
path.append(resultsname+'/' + casename + '/field_plot_data' )
path.append(resultsname+'/' + casename + '/line_plot_data' )
for num,path_id in enumerate(path):
    isExists = os.path.exists(path_id)
    if not isExists:
        os.mkdir(path_id)   # Make a directory for the results.

if print_convergence == 'off':
    set_log_level(40)    # Suppress Newton convergence information if desired.
print('casename is ' + casename)

# ======================= Create VTK files ==========================
vtkfile_u = File(resultsname+'/'+casename+'/field_plot_data/u.pvd')            # velocity
vtkfile_p = File(resultsname+'/'+casename+'/field_plot_data/p.pvd')            # pressure
vtkfile_sigma = File(resultsname+'/'+casename+'/field_plot_data/sigma.pvd')    # stress
vtkfile_strain_rate = File(resultsname+'/'+casename+'/field_plot_data/strain_rate.pvd')        # strain
fname = open('residual.txt','w')

# ======================= import the mesh ==========================
# import the xml mesh
mesh = Mesh('./meshes/'+'tides'+'_DX'+str(int(DX_s))+\
            '_Lngth'+str(int(Lngth))+'_Slope0_02'+'_U'+"{:0>2d}".format(int(U0*3.154e7))+'.xml')
new_mesh = File(resultsname+'/'+casename+'/tides'+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth)) + '_Slope0_02'\
            +'_U'+"{:0>2d}".format(int(U0*3.154e7))+'.xml')

# Define arrays for saving surfaces, lake volume, water pressure, and
# grounding line positions over time.
Gamma_s = np.zeros((nx,nt))       # Basal surface
Gamma_h = np.zeros((nx,nt))       # Upper surface
s_mean = np.zeros(nt)             # Mean elevation of ice-water interface
h_mean = np.zeros(nt)             # Mean elevation of surface above ice-water interface
x_left = np.zeros(nt)             # Left grounding line position
x_right = np.zeros(nt)            # Right grounding line position
P_res = np.zeros(nt)              # Penalty functional residual
t = 0                             # Time

# ======================= Begin time stepping ==========================
for i in range(nt):
    print('-----------------------------------------------')
    print('Iteration '+str(i+1)+' out of '+str(nt))

    # Set initial conditions.
    if t==0:   
        s_mean_i = s_mean0                    # Mean ice-water elevation.
        F_h = lambda x: Hght                  # Ice-air surface function
        F_s = lambda x: interface(x)          # Lower surface function

        if model == 'marine':
            w_0 = get_zero_m(mesh)              # Placeholder for first iteration,
                                              # used for computing surface elevation functions
        mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w_0,mesh,dt,interface,surface)

    # Solve the Stokes problem.
    # Returns solutions "w" and penalty functional residual "P_res_i"
    if t==0 and model == 'marine' and tides == 'on':
        fFile = HDF5File(MPI.comm_world,"w_init_DX"+str(int(DX_s))+"_L"+str(int(Lngth))+"_Slope0_02"\
        +'_U'+"{:0>2d}".format(int(U0*3.154e7))+'.h5',"r")
        fFile.read(w_0,"/f")
        fFile.close()

    # Deep copy of the mesh at last step
    w,P_res_i,_strain_rate = stokes_solve_marine(mesh,F_h,t,w_0)
    w.set_allow_extrapolation(True)

    # save the old mesh and solution for next step calculation
    mesh_0 = Mesh(mesh)
    w_0 = get_zero_m(mesh_0)
    w_0.vector().set_local(w.vector().get_local())
    w_0.set_allow_extrapolation(True)

    # Solve the surface kinematic equations, move the mesh, and compute the
    # grounding line positions.

    mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w,mesh,dt,F_s,F_h)

    # Save quantities of interest.
    P_res[i] = P_res_i                        # residual                      
    s_mean[i] = s_mean_i                      
    h_mean[i] = h_mean_i
    x_left[i] = XL                            # left grounding line
    x_right[i] = XR                           # right grounding line
    Gamma_s[:,i] = F_s(X_fine)                # bottom boundary
    Gamma_h[:,i] = F_h(X_fine)                # surface boundary

    # Save Stokes solution
    _u, _p, _sigma = w.split()
    _u.rename("vel", "U")
    _p.rename("press","P")
    _sigma.rename("stress","SIGMA")
    _strain_rate.rename("strain","STRAIN")

    vtkfile_u << (_u,t)
    vtkfile_p << (_p,t)
    vtkfile_sigma << (_sigma,t)
    vtkfile_strain_rate << (_strain_rate,t)
    # Vtkfile_stress_rate << (_stress_rate,t)

    # Update time
    t += dt
    # Save the checkpoint
    if model == 'marine':
        fFile = HDF5File(MPI.comm_world,resultsname+'/'+casename+"/w_init_DX"+str(int(DX_s))+"_L"+str(int(Lngth))+".h5","w")
        fFile.write(w_0,"/f")
        fFile.close()

    # Print information of interest.
    print('Left grounding line: '+str(x_left[i]/1000.0)+' km')
    print('Right grounding line: '+str(x_right[i]/1000.0)+' km')

# Save quantities of interest.
t_arr = np.linspace(0,t_final,num=int(nt_per_year*t_final/3.154e7))

# save the mesh when initiating the mesh
if model == 'marine' and tides == 'on':
    new_mesh << mesh

np.savetxt(resultsname+'/'+casename+'/line_plot_data/Gamma_s',Gamma_s)    # see definition above
np.savetxt(resultsname+'/'+casename+'/line_plot_data/Gamma_h',Gamma_h)    # "   "
np.savetxt(resultsname+'/'+casename+'/line_plot_data/s_mean',s_mean)      # ...
np.savetxt(resultsname+'/'+casename+'/line_plot_data/h_mean',h_mean)
np.savetxt(resultsname+'/'+casename+'/line_plot_data/x_left',x_left)
np.savetxt(resultsname+'/'+casename+'/line_plot_data/x_right',x_right)    # ...
np.savetxt(resultsname+'/'+casename+'/line_plot_data/P_res',P_res)        # "   "
np.savetxt(resultsname+'/'+casename+'/line_plot_data/X',X_fine)           # X = spatial coordinate
np.savetxt(resultsname+'/'+casename+'/line_plot_data/t',t_arr)            # t = time coordinate