#-------------------------------------------------------------------------------
# This file contains functions that:
# (1) define the boundaries (ice-air,ice-water,ice-bed) of the mesh,
# (2) mark the boundaries of the mesh, AND ...
# (3) create Dirichlet boundary conditions on one or both side walls of the domain.
#-------------------------------------------------------------------------------
from params import tol,U0,model,Lngth,Hght,bed_slope,A0,rho_i,g,bed_slope,n
from geometry import bed
import numpy as np
from dolfin import *

#-------------------------------------------------------------------------------
# Define SubDomains for ice-water boundary, ice-bed boundary, inflow (x=0) and
# outflow (x=Length of domain). The parameter 'tol' is a minimal water depth
# used to distinguish the ice-water and ice-bed surfaces.

class WaterBoundary(SubDomain):
    # Ice-water boundary.
    # Note: This boundary is marked first and all of the irrelevant portions are
    # overwritten by the other boundary markers. This results in a "last grounded"
    # scheme as described in Gagliardini et al. (2016), The Cryosphere.
    def inside(self, x, on_boundary):
        return (on_boundary and (x[1]<0.5*Hght))

class BedBoundary(SubDomain):
    # Ice-bed boundary
    def inside(self, x, on_boundary):
        return (on_boundary and ((x[1]-bed(x[0]))<=tol))

class LeftBoundary(SubDomain):
    # Left boundary
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[0])<tol)

class RightBoundary(SubDomain):
    # Right boundary
    def inside(self, x, on_boundary):
        return (on_boundary and np.abs(x[0]-Lngth)<tol)

#-------------------------------------------------------------------------------

def mark_boundary(mesh):
    # Assign markers to each boundary segment (except the upper surface).
    # This is used at each time step to update the markers.
    #
    # Boundary marker numbering convention:
    # 1 - Left boundary
    # 2 - Right boundary
    # 3 - Ice-bed boundary
    # 4 - Ice-water boundary
    # 5 - Ice-supraglacial boundary

    boundary_markers = MeshFunction('size_t', mesh,dim=1)
    boundary_markers.set_all(0)

    # Mark ice-water boundary
    bdryWater = WaterBoundary()
    bdryWater.mark(boundary_markers, 4)

    # Mark ice-bed boundary
    bdryBed = BedBoundary()
    bdryBed.mark(boundary_markers, 3)

    # Mark inflow boundary
    bdryLeft = LeftBoundary()
    bdryLeft.mark(boundary_markers, 1)

    # Mark outflow boundary
    bdryRight = RightBoundary()
    bdryRight.mark(boundary_markers, 2)

    return boundary_markers

#------------------------------------------------------------------------------

def apply_bcs(W,Fh,boundary_markers):
    # Apply inflow and outflow boundary conditions to the system.
    # These are applied to the horizontal velocity component.
    angle = np.arctan(bed_slope)
    # linear viscosity

    # Dirichlet BC: inflow velocity set-up
    # u_left = Expression('+scos*(u_surf-2.0*A/(n+1)*(pow(rho_i*g*ssin,n))*pow((h-x[1])*scos,n+1))',\
    #     degree=2, scos=np.cos(angle), ssin=np.sin(angle), n=int(n), u_surf=U0, A=A0,rho_i=rho_i, g=g, h=float(Fh(0)))
    # v_left = Expression('-ssin*(u_surf-2.0*A/(n+1)*(pow(rho_i*g*ssin,n))*pow((h-x[1])*scos,n+1))',\
    #     degree=2, scos=np.cos(angle), ssin=np.sin(angle), n=int(n), u_surf=U0, A=A0,rho_i=rho_i, g=g, h=float(Fh(0)))

    bcu1 = DirichletBC(W.sub(0).sub(0), Constant(U0), boundary_markers,1)
    # bcu1 = DirichletBC(W.sub(0).sub(0), u_left, boundary_markers,1)
    bcu2 = DirichletBC(W.sub(0).sub(0), Constant(U0), boundary_markers,2)
    # bcu3 = DirichletBC(W.sub(0).sub(1), v_left, boundary_markers,1)

    if model == 'lake':
        BC = [bcu1,bcu2]
    elif model == 'marine':
        BC = [bcu1]
    return BC
