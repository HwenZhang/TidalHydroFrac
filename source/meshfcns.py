#------------------------------------------------------------------------------
# These functions are used to:
# (1) update the mesh at each timestep by solving the
#     surface kinematic equations, AND...
# (2) compute the grounding line positions.
#------------------------------------------------------------------------------

from params import tol,Lngth,Hght,X_fine
from dolfin import *
from fenics import *
import numpy as np
from geometry import bed
from scipy.interpolate import interp1d

# ------------------------------------------------------------------------------

def mesh_routine(w,mesh,dt,F_s,F_h):
    # This function solves the surface kinematic equations and moves the mesh.
    # The mean elevation of the ice-water and ice-air interfaces are also
    # computed and returned.

    # Compute slopes of free surfaces.
    # Returns: (1) FEniCS functions baseslope and surfslope, AND...
    #          (2) Python functions (F_h, F_s) of the surface elevations.
    
    baseslope,F_s = get_baseslope(mesh,F_s)
    surfslope, F_h = get_surfslope(mesh,F_h)

    # Get maximum and minimum grounding line positions
    XL,XR = get_glines(F_s)

    # Move the mesh
    move_mesh(mesh,baseslope,surfslope,dt,F_s,F_h,w)

    # Compute mean elevation of ice-water and ice-air surfaces.
    s_mean = np.mean(F_s(X_fine)[F_s(X_fine)-tol>bed(X_fine)])
    h_mean = np.mean(F_h(X_fine)[F_s(X_fine)-tol>bed(X_fine)])

    return mesh,F_s,F_h,s_mean,h_mean,XL,XR


#------------------------------------------------------------------------------
def move_mesh(mesh,baseslope,surfslope,dt,F_s,F_h,w):
    # This function computes the surface displacements and moves the mesh.

    # print(w.get_allow_extrapolation())
    M = mesh.coordinates()                           # Get mesh coordinates.

    w0 = w.sub(0).sub(1).compute_vertex_values(mesh) # Get vertical velocity at nodes.
    u0 = w.sub(0).sub(0).compute_vertex_values(mesh) # Get horizontal velocity at nodes.

    sx = baseslope.compute_vertex_values(mesh)       # Get lower surface slope at nodes.
    hx = surfslope.compute_vertex_values(mesh)       # Get upper surface slope at nodes.

    # Compute vertical displacements via the kinematic equation:
    # dZ/dt = w - u * dZ/dx
    # for Z = s(x,t) and Z = h(x,t).
    disp0 = w0 - u0*sx                               # Compute lower surface displacement.
    disp1 = w0 - u0*hx                               # Compute upper surface displacement.

    # Mark all of the vertices on the boundary of the mesh
    V = FunctionSpace(mesh, 'CG', 1)
    bc = DirichletBC(V, 1, DomainBoundary())
    u = Function(V)
    bc.apply(u.vector())
    d2v = dof_to_vertex_map(V)
    vertices_on_boundary = d2v[u.vector() == 1]
    vertices_on_boundary = np.sort(vertices_on_boundary)

    # Loop over nodes in the boundary and displace them vertically according to
    # the velocity solution and surface slope.
    for i in vertices_on_boundary:
        
        # BOTTOM surface: ice-water interface
        if np.abs(M[i,1]-F_s(M[i,0]))<tol:

            M[i,1] += dt*disp0[i]
            # print(dt*disp0[i])
            # If new y-value is below the bed, set equal to the bed elevation
            if M[i,1]-bed(M[i,0])<0:
                M[i,1] = bed(M[i,0])

        #Top surface: ice-air interface
        elif np.abs(M[i,1]-F_h(M[i,0]))<tol:
            M[i,1] += dt*disp1[i]

#------------------------------------------------------------------------------
def get_baseslope(mesh,F_s):
    # This function computes the slope of the lower surface and returns it as a FEniCS function.

    # Get x and y values of boundary nodes
    bmesh = BoundaryMesh(mesh,'exterior')
    M = bmesh.coordinates()

    # (M[:,0]>tol) & (M[:,0]<Lngth-tol) : exclude both sides
    # abs(M[:,1]-interface(M[:,0]))<tol : exclude the top
    X =  M[:,0][ (M[:,0]>tol) & (M[:,0]<Lngth-tol) & (abs(M[:,1]-F_s(M[:,0]))<0.1*Hght)]
    Y =  M[:,1][ (M[:,0]>tol) & (M[:,0]<Lngth-tol) & (abs(M[:,1]-F_s(M[:,0]))<0.1*Hght)]

    Y = Y[np.argsort(X)]
    X = np.sort(X)

    s_left = np.min(M[:,1][M[:,0]<tol])
    s_right = np.min(M[:,1][np.abs(M[:,0]-Lngth)<tol])

    # Append values at x=0 and x=Lngth.
    Y = np.append(Y,s_right)
    Y = np.insert(Y,0,s_left)

    X = np.append(X,Lngth)
    X = np.insert(X,0,0)

    # Use SciPy to interpolate the lower surface
    F_s = interp1d(X,Y,kind='linear',fill_value='extrapolate',bounds_error=False)

    # Define a FEniCS expression for the lower surface elevation
    class ExpressionPhi_s(UserExpression):
        def eval(self,value,x):
            value[0] = F_s(x[0])

    # Compute the slope of the lower surface in FEniCS
    V = FunctionSpace(mesh,'CG',1)
    Phi_s = ExpressionPhi_s(element=V.ufl_element(),domain=mesh)

    sx = Dx(Phi_s,0)
    baseslope = Function(V)
    baseslope.assign(project(sx,V))

    return baseslope, F_s

#------------------------------------------------------------------------------

def get_surfslope(mesh,F_h):
    # This function computes the upper surface slope

    # Get coordinates of mesh nodes on boundary.
    bmesh = BoundaryMesh(mesh,'exterior')
    M = bmesh.coordinates()

    X =  M[:,0][ (M[:,0]>tol) & (M[:,0]<Lngth-tol) & (abs(M[:,1]-F_h(M[:,0]))<0.1*Hght)]
    Y =  M[:,1][ (M[:,0]>tol) & (M[:,0]<Lngth-tol) & (abs(M[:,1]-F_h(M[:,0]))<0.1*Hght)]

    Y = Y[np.argsort(X)]
    X = np.sort(X)

    h_left = np.max(M[:,1][M[:,0]<tol])
    h_right = np.max(M[:,1][np.abs(M[:,0]-Lngth)<tol])

    # Append values at x=0 and x=Lngth.
    Y = np.append(Y,h_right)
    Y = np.insert(Y,0,h_left)

    X = np.append(X,Lngth)
    X = np.insert(X,0,0)

    # Interpolate the boundary points:
    F_h = interp1d(X,Y,kind='linear',fill_value='extrapolate',bounds_error=False)

    # Define a FEniCS expression for the upper surface elevation
    class ExpressionPhi_h(UserExpression):
        def eval(self,value,x):
            value[0] = F_h(x[0])

    # Compute slope of upper surface
    V = FunctionSpace(mesh,'CG',1)
    Phi_h = ExpressionPhi_h(element=V.ufl_element(),domain=mesh)

    hx = Dx(Phi_h,0)
    surfslope = Function(V)
    surfslope.assign(project(hx,V))

    return surfslope,F_h

#------------------------------------------------------------------------------

def get_glines(F_s):
    # Computes minimum and maximum grounding line positions given the
    # lower surface elevation function s and the bed geometry.

    s = F_s(X_fine)
    s_new = np.copy(s)

    x = np.copy(X_fine)
    x_new = np.copy(X_fine)

    key = 1e10

    # Mark points on ice-bed boundary (by assigning a ridiculously large value： 1e10)
    s[s-bed(x)<tol] = key
    s_new[s_new-bed(x)<tol] = key

    # Loop over points on lower boundary
    for j in range(np.shape(s)[0]-1):
        if s[j+1] < 0.9*key and  s[j-1] < 0.9*key:
        # If *both* neighboring points are on the ice-water boundary, then
        # mark these points! (these cannot be grounding lines)
            s_new[j] = key

    # Mark last point
    if s[-2] < 0.9*key:
        s_new[-1] = key

    # All points except grounding lines have now been marked.
    glines = x_new[s_new<0.9*key]
    XL = np.min(glines) # Minimum grounding line
    XR = np.max(glines) # Maximum grounding line

    return XL,XR