#-------------------------------------------------------------------------------
# Generate 2D domain for marine ice sheet problem.
# Modified by Hanwen Zhang on 1st Aug 2022.

# The mesh is generated as a .xdmf file based on gmsh. 
#-------------------------------------------------------------------------------
import sys
import os
sys.path.insert(0, './source')
sys.path.insert(0, './scripts')
import gmsh
import meshio
import numpy as np
from geometry import interface
from params import Hght,Lngth,model,DX_s,DX_h,bed_slope

nx = 6                                  # Number of grid points in x direction.
X = np.linspace(0,Lngth,num=nx)         # array for horizontal coordinate lower
S = interface(X)                        # Ice-water interface array
alpha = np.arctan(bed_slope)            # bed slope angle

gmsh.initialize(sys.argv)
gmsh.model.add("mesh")
opengmsh=True

fname = model+'_DX'+str(int(DX_s))+'_Lngth'+str(int(Lngth))+'_Slope2e_2'

if not os.path.exists('./meshes/' + fname):
    os.mkdir('./meshes/' + fname)   # Make a directory for the results.
#--------------------- points ------------------------ 
point=[]
# Bottom interface
# index from 0 to nx-1
point.append(gmsh.model.geo.addPoint(X[0],S[0],0,DX_h/2,0))
point.append(gmsh.model.geo.addPoint(Lngth/2,interface(Lngth/2),0,DX_s,1))
point.append(gmsh.model.geo.addPoint(X[nx-1],S[nx-1],0,DX_h/2,2))

point.append(gmsh.model.geo.addPoint(X[nx-1],S[nx-1]+Hght,0,DX_h,3))
point.append(gmsh.model.geo.addPoint(Lngth/2+Hght*(1/np.sin(alpha)-1/np.tan(alpha)),Hght,0,DX_h,4))
point.append(gmsh.model.geo.addPoint(X[0],S[0]+Hght/np.cos(alpha),0,DX_h,5))


#--------------------- lines ------------------------ 
lines = []
for i in range(nx-1):
    lines.append(gmsh.model.geo.addLine(i, i+1, i))
lines.append(gmsh.model.geo.addLine(nx-1, 0, nx-1))
loop = list(range(nx))

# boundary loop and 2d plane
boundary_loop = gmsh.model.geo.addCurveLoop(loop,nx)
plane_surface = gmsh.model.geo.addPlaneSurface([nx], 1)
gmsh.model.geo.synchronize()

# boundary markers
boundary0 = gmsh.model.addPhysicalGroup(1, [0], 3)
gmsh.model.setPhysicalName(1, boundary0, "BedBdry")

boundary1 = gmsh.model.addPhysicalGroup(1, [1], 4)
gmsh.model.setPhysicalName(1, boundary1, "WaterBdry")

boundary2 = gmsh.model.addPhysicalGroup(1, [2], 2)
gmsh.model.setPhysicalName(1, boundary2, "OutflowBdry")

boundary3 = gmsh.model.addPhysicalGroup(1, [3,4], 0)
gmsh.model.setPhysicalName(1, boundary3, "TopBdry")

boundary4 = gmsh.model.addPhysicalGroup(1, [5], 1)
gmsh.model.setPhysicalName(1, boundary4, "InflowBdry")

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
meshio.write('./meshes/'+fname+'/'+fname+'.xdmf',triangle_mesh)


line_mesh = meshio.Mesh(points=mesh_from_file.points,
                      cells=[("line", mesh_from_file.get_cells_type("line"))],
                      cell_data={"Boundaries": [mesh_from_file.cell_data_dict['gmsh:geometrical']['line']]},
                      field_data=mesh_from_file.field_data)
line_mesh.prune_z_0()
meshio.write('./meshes/'+fname+'/'+fname+'_facet.xdmf',line_mesh)