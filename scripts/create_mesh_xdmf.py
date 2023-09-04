#-------------------------------------------------------------------------------
# Generate 2D domain for subglacial lake problem.
# Modified by Hanwen Zhang on 1st Aug 2022.
# Generate .xdmf file based on gmsh

# Note: This always assumes a flat upper surface.
#-------------------------------------------------------------------------------
import sys
import os
sys.path.insert(0, './source')
import gmsh
import meshio
import numpy as np
import matplotlib.pyplot as plt
from geometry import bed,interface
from params import Hght,Lngth,model,DX_s,DX_h,model_setup
import subprocess

nx = int(Lngth/DX_s)+1                  # Number of grid points in x direction.
nxh = int(Lngth/DX_h)+1                  # Number of grid points in x direction.
X = np.linspace(0,Lngth,num=nx)         # array for horizontal coordinate lower
Xh = np.linspace(0,Lngth,num=nxh)         # array for horizontal coordinate upper
S = interface(X)                        # Ice-water interface array
H = interface(Xh)                        # Ice-water interface array

gmsh.initialize(sys.argv)
gmsh.model.add("mesh")
opengmsh=True

fname = model+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope0_05'
#--------------------- points ------------------------ 
point=[]
# Bottom interface
# index from 0 to nx-1
for i in range(nx):
    point.append(gmsh.model.geo.addPoint(X[i],S[i],0,DX_s,i))

# Upper surface
for i in range(nxh):
    point.append(gmsh.model.geo.addPoint(Xh[nxh-1-i],H[nxh-1-i]+Hght,0,DX_h,nx+i))

#--------------------- lines ------------------------ 
lines = []
for i in range(nxh+nx-1):
    lines.append(gmsh.model.geo.addLine(i, i+1, i))
lines.append(gmsh.model.geo.addLine(nxh+nx-1, 0, nxh+nx-1))
loop = list(range(nxh+nx))
# boundary loop and 2d plane

boundary_loop = gmsh.model.geo.addCurveLoop(loop,nxh+nx)
plane_surface = gmsh.model.geo.addPlaneSurface([nxh+nx], 1)
gmsh.model.geo.synchronize()

# Mesh Size Initiation
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.model.mesh.generate(2)
gmsh.write(fname+'.msh2')

# change file name
os.rename(fname+'.msh2', fname+'.msh')
if opengmsh:
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
gmsh.finalize()

#============================================================
# convert msh file to xdmf file
mesh_from_file = meshio.read(fname+'.msh')
prune_z=True
triangle_mesh = meshio.Mesh(points=mesh_from_file.points,
                      cells=[("triangle", mesh_from_file.get_cells_type("triangle"))],
                      cell_data={"Subdomain": [mesh_from_file.cell_data_dict['gmsh:physical']['triangle']]},
                      field_data=mesh_from_file.field_data)
triangle_mesh.prune_z_0()
meshio.write(fname+'_mesh.xdmf',triangle_mesh)

line_mesh = meshio.Mesh(points=mesh_from_file.points,
                      cells=[("line", mesh_from_file.get_cells_type("line"))],
                      cell_data={"Boundaries": [mesh_from_file.cell_data_dict['gmsh:geometrical']['line']]},
                      field_data=mesh_from_file.field_data)
line_mesh.prune_z_0()
meshio.write(fname+'_facet_mesh.xdmf',line_mesh)

#------------------Generate initial ice-water/ice-bed interface-----------------
def top_surface(x):
    if model == 'lake':
        Int = 0.5*(bed(x) - 5 + np.abs(bed(x) - (-5)))

        if model_setup == 'wedge_test':
            Int = 0.5*(bed(x) - 7.5 + np.abs(bed(x) - (-7.5)))


    elif model == 'marine':
        Int = 0.5*(bed(x)  + np.abs(bed(x) ))
    return Int
#-------------------------------------------------------------------------------

