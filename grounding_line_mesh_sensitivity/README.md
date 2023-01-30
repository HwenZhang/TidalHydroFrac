# TidalHydroFrac

FEniCS scipts to model the tidal response of the Amery Ice Shelf. Below is the schematic showing the model set-up.
![Image text](https://github.com/HwenZhang/TidalHydroFrac/blob/147148f5916b7197c94a07abe23951a49d448c2f/grounding_line_mesh_sensitivity/image/schematic.png)


# Ice Rheology
## Incompressible Viscoelasticity
  nonlinear viscosity $\eta$: A=3e-24, n=3
  shear modulus: $\mu=0.18\times 10^{9}$ $\text{Pa}$

# Boundary condition
## Left boundary 
Constant inflow: $U_0 = 9$ m/yr in x-direction.
## Right boundary
Ice overburden stress
$$\sigma_{xx}=\rho_i g \left(H-z\right)$$

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
