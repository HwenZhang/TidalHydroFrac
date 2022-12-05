#-------------------------------------------------------------------------------
# This program solves a nonlinear Stokes flow problem describing groundning line
# dynamics in two settings:
#
# (1) A marine ice sheet, OR ...
# (2) A subglacial lake that is (potentially) undergoing water volume changes.
#
# Choose model setup (1) or (2) in the params.py file.
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
                    print_convergence,X_fine,nx,tides,DX_s,model_setup)
from params import (rho_i,g,tol,B,rm2,rho_w,C,eps_p,eps_v,sea_level,dt,slope_str,
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

# ======================= Create VTK files ==========================
vtkfile_u = File(resultsname+'/'+casename+'/field_plot_data/u.pvd')            # velocity
vtkfile_p = File(resultsname+'/'+casename+'/field_plot_data/p.pvd')            # pressure
vtkfile_sigma = File(resultsname+'/'+casename+'/field_plot_data/sigma.pvd')    # stress
vtkfile_eta = File(resultsname+'/'+casename+'/field_plot_data/eta.pvd')        # velocity
vtkfile_res = File(resultsname+'/'+casename+'/field_plot_data/stress_res.pvd')        # residual
fname = open('residual.txt','w')


# Load and read mesh
if model == 'marine' and tides=='on':
    # tidal simulation
    meshname = 'tides'+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope'+slope_str+'.xdmf'
    # import the mesh
    mesh = Mesh('tides_DX50_Lngth20000_Slope0_05.xml')

elif model == 'marine' and tides=='off':
    # mesh initiation
    meshname = model+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope'+slope_str+'.xdmf'
    meshdict = model+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope'+slope_str
    new_mesh = File('tides'+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope'+slope_str+'.xml')

    # import the mesh
    mesh = Mesh()
    with XDMFFile(MPI.comm_world, "./meshes/"+meshdict+'/'+meshname) as meshfile:
        meshfile.read(mesh)
        mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    print("The MeshValueCollection: ", mvc)
else:
    exit()

# Define arrays for saving surfaces, lake volume, water pressure, and
# grounding line positions over time.
Gamma_s = np.zeros((nx,nt))       # Basal surface
Gamma_h = np.zeros((nx,nt))       # Upper surface
s_mean = np.zeros(nt)             # Mean elevation of ice-water interface
h_mean = np.zeros(nt)             # Mean elevation of surface above ice-water interface
x_left = np.zeros(nt)             # Left grounding line position
x_right = np.zeros(nt)            # Right grounding line position
P_res = np.zeros(nt)              # Penalty functional residual
t = 0                             # time

# Begin time stepping
for i in range(nt):
    print('-----------------------------------------------')
    print('Iteration '+str(i+1)+' out of '+str(nt))

    if t==0:
        # Set initial conditions.
        s_mean_i = s_mean0                    # Mean ice-water elevation.
        F_h = lambda x: Hght                  # Ice-air surface function
        F_s = lambda x: interface(x)          # Lower surface function
 
        if model == 'marine':
            w = get_zero_m(mesh)              # Placeholder for first iteration,
                                              # used for computing surface elevation functions
        mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w,mesh,dt,interface,surface)

    # Solve the Stokes problem.
    # load the initial guess saved before
    # if t==0:
    #     fFile = HDF5File(MPI.comm_world,"w_init_DX"+str(int(DX_s))+"_L"+str(int(Lngth))+".h5","r")
    #     fFile.read(w,"/f")
    #     fFile.close()
    # Returns solutions "w" and penalty functional residual "P_res_i"
    w,P_res_i,_eta = stokes_solve_marine(mesh,F_h,t,w)

    # Solve the surface kinematic equations, move the mesh, and compute the
    # grounding line positions.
    mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w,mesh,dt,F_s,F_h)
    # Save quantities of interest.
    P_res[i] = P_res_i
    s_mean[i] = s_mean_i
    h_mean[i] = h_mean_i
    x_left[i] = XL
    x_right[i] = XR
    Gamma_s[:,i] = F_s(X_fine)
    Gamma_h[:,i] = F_h(X_fine)

    # Save Stokes solution
    _u, _p, _sigma = w.split()
    _u.rename("vel", "U")
    _p.rename("press","P")
    _sigma.rename("stress","SIGMA")
    _eta.rename("viscosity","ETA")

    vtkfile_u << (_u,t)
    vtkfile_p << (_p,t)
    vtkfile_sigma << (_sigma,t)
    vtkfile_eta << (_eta,t)
    
    # Update time
    t += dt

    # Print information of interest.
    print('Left grounding line: '+str(x_left[i]/1000.0)+' km')
    print('Right grounding line: '+str(x_right[i]/1000.0)+' km')

# Save quantities of interest.
t_arr = np.linspace(0,t_final,num=int(nt_per_year*t_final/3.154e7))

if model == 'marine' and tides == 'off':
    XDMFFile(MPI.comm_world, "./meshes/"+'tides'+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope'+slope_str+'.xdmf').write(mesh)
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