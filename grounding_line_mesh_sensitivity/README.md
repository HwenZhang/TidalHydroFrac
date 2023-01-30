# TidalHydroFrac

FEniCS scipts to model the tidal response of the Amery Ice Shelf.

# Ice Rheology
## Incompressible Viscoelasticity
  nonlinear viscosity $\eta$: A=3e-24, n=3
  shear modulus: $\mu=0.18\times 10^{9}$ $\text{Pa}$

# Boundary condition
## Left boundary 
Constant inflow: $U_0 = 9$ m/yr in x-direction.
## Right boundary
Ice overburden stress
\begin{equation}
  \sigma_{xx}=\rho_i g \left(H-z\right)
\end{equation}

## Top boundary
Free surface.

## Bottom boundary
### Ice-bedrock interface
Penalty method.

### Ice-ocean interface
Hydrostatic pressure.
  
# Geometry
Length: $L=20$ km
Height: $H=500$ m
Tidal amplitude: $A=1.0$ m