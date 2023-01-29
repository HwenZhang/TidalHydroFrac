# This file contains the functions needed for solving the Stokes system.
from params import rho_i,g,tol,B,rm2_s,rm2,rho_w,C,eps_p,eps_v,sea_level,dt,quad_degree,Lngth,U0,mu,eps_vs
from boundaryconds import mark_boundary,apply_bcs
from hydrology import Vdot, sl_change
import numpy as np
from dolfin import *
from fenics import *

def dPi(u,nu):
        # Derivative of penalty functional for enforcing impenetrability on the ice-bed boundary.
        un = dot(u,nu)
        return un+abs(un)

def Pi(u,nu):
        # Penalty functional for enforcing impenetrability on the ice-bed boundary.
        un = dot(u,nu)
        return 0.5*(un**2.0+un*abs(un))

def weak_form_marine(u,p,tau,u_old,p_old,tau_old,u_test,p_test,tau_test,f,g_base,g_in,g_out,ds,nu,T):
    # Weak form of the marine ice sheet problem
    # viscosity (non-newtonian fluid)
    eta = B*(inner(sym(grad(u)),sym(grad(u)))+Constant(eps_v))**(rm2/2.0)

    # strain rate
    epsilon = 0.5*(nabla_grad(u)+nabla_grad(u).T)
    # SUPG stabilisation unused
    # relaxation time
    lamda = eta / mu
    # incompressibiity
    w1 = p_test*div(u)*dx
    # stokes equation
    w2 = (-div(u_test)*p+inner(tau,grad(u_test)))*dx - inner(f, u_test)*dx\
         + (g_base+Constant(rho_w*g*dt)*inner(u,nu))*inner(nu, u_test)*(ds(3)+ds(4))\
         + Constant(C)*((inner(dot(T,u),dot(T,u))+Constant(eps_v))**(rm2_s/2.0))*inner(dot(T,u),dot(T,u_test))*ds(3)\
         + g_out*dot(nu, u_test)*ds(2)\
         + Constant(1.0/eps_p)*dPi(u,nu)*dot(u_test,nu)*ds(3)
      
    # constitutive relation (upper-convected rate)     
    w3 = inner(tau+lamda*((tau-tau_old)/dt+dot(u,nabla_grad(tau))-dot(grad(u),tau)-dot(tau, grad(u).T))\
         -2.*eta*epsilon,tau_test)* dx   
    Fw = w1 + w2 + w3
    return Fw


def stokes_solve_marine(mesh,F_h,t,w_old):
        # Stokes solver using Taylor-Hood elements.

        # Define function spaces
        P1 = FiniteElement('CG',triangle,1)     # Pressure
        P2 = VectorElement('CG',triangle,2)     # Velocity
        P3 = TensorElement("CG",triangle,2)    # Stress

        element = MixedElement([P2,P1,P3])
        W = FunctionSpace(mesh,element)        # Function space for (u,p)
        # Assign the function w from previous mesh to the current mesh
        (u_old,p_old,tau_old) = w_old.split(True)
        #---------------------Define variational problem------------------------
        # project the velocity field to the new mesh
        w_prev = project(w_old,W) # used for compute \pdiff{\tau}{t}
        # w_prev = interpolate(w_old,W) # used for compute \pdiff{\tau}{t}
        (u_prev,p_prev,tau_prev) = w_prev.split(True)
        # copy the previous step solution w_old to w as an initial guess
        w = Function(W)
        dw = TrialFunction(W)
        w.vector().set_local(w_old.vector().get_local())
        (u,p,tau) = split(w)

        # TestFunctions
        (u_test,p_test,tau_test) = TestFunctions(W)
        # Neumann condition at outflow boundary
        h_out = float(F_h(Lngth))        # Surface elevation at outflow boundary
        h_in = float(F_h(0))        # Surface elevation at outflow boundary
        g_out = Expression('rho_i*g*(h_out-x[1])',rho_i=rho_i,g=g,h_out=h_out,degree=1)
        g_in = Expression('rho_i*g*(h_in-x[1])',rho_i=rho_i,g=g,h_in=h_in,degree=1)
        # Neumann condition at ice-water boundary
        g_base = Expression('rho_w*g*(sea_level-x[1])',rho_w=rho_w,g=g,sea_level=sea_level+sl_change(t),degree=1)

        f = Constant((0,-rho_i*g))        # Body force
        nu = FacetNormal(mesh)            # Outward-pointing unit normal to the boundary
        I = Identity(2)                   # Identity tensor
        T = I - outer(nu,nu)              # Tangential projection operator

        # Mark bounadries of mesh and define a measure for integration
        boundary_markers= mark_boundary(mesh)
        ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

        # Define weak form and apply boundary conditions on the inflow boundary
        bcs_u =  apply_bcs(W,F_h,boundary_markers)    # Apply Dirichlet BC

        # element size
        # max_h = CellDiameter(mesh)
        # Solve for (u,p).
        Fw = weak_form_marine(u,p,tau,u_prev,p_prev,tau_prev,u_test,p_test,tau_test,f,g_base,g_in,g_out,ds,nu,T)

        #solve(Fw == 0, w, bcs=bcs_u,solver_parameters={"newton_solver":{"relative_tolerance": 1e-14,"absolute_tolerance": 5e-4,\
        #        "maximum_iterations":100,"relaxation_parameter":1.0}},\
        #                form_compiler_parameters={"quadrature_degree":quad_degree,"optimize":True,"eliminate_zeros":False})
        parameters["form_compiler"]["quadrature_degree"]=quad_degree
        parameters["form_compiler"]["optimize"]=True

        jacobian = derivative(Fw, w, dw)
        pde_nonlinear = NonlinearVariationalProblem(Fw, w, bcs_u,jacobian)
        solver_nonlinear = NonlinearVariationalSolver(pde_nonlinear)
        prm = solver_nonlinear.parameters
        newton_prm = prm['newton_solver']

        newton_prm['linear_solver'] = 'umfpack'
        newton_prm['krylov_solver']['maximum_iterations'] = 2000
        newton_prm['krylov_solver']['nonzero_initial_guess'] = True
        newton_prm['krylov_solver']['monitor_convergence'] = False

        # relaxation parameter
        newton_prm['relaxation_parameter'] = 0.40

        newton_prm["relative_tolerance"] = 1e-14
        newton_prm["absolute_tolerance"] = 5e-3
        newton_prm["maximum_iterations"] = 150
        newton_prm['error_on_nonconvergence'] = True   #########!!!!!
        solver_nonlinear.solve()

        # Compute penalty functional residiual
        P_res = assemble(Pi(u,nu)*ds(3))
        # Save the viscosity
        eta_space = FunctionSpace(mesh,P1)
        eta = project(B*(inner(sym(grad(u)),sym(grad(u)))+Constant(eps_v))**(rm2/2.0),eta_space)
        return w,P_res,eta

def get_zero_m(mesh):
        # Get zero element of function space for marine ice sheet problem.
        # Only used for setting initial conditions; see main.py.

        # Define function spaces
        P1 = FiniteElement('CG',triangle,1)     # Pressure
        P2 = VectorElement('CG',triangle,2)     # Velocity
        P3 = TensorElement("CG",triangle,2)     # Stress
        element = MixedElement([P2,P1,P3])
        W = FunctionSpace(mesh,element)         # Function space for (u,p)
        w = Function(W)
        return w