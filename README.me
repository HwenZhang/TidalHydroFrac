# TidalHydroFrac

FEniCS scipts to model the tidal response of the Amery Ice Shelf. Below is the schematic showing the model set-up. The code is modified from the variational formulation proposed by Stubblefield et al., 2019. Below attached the link to the paper:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/variational-formulation-of-marine-icesheet-and-subglaciallake-groundingline-dynamics/70AA4C7F565E3C8B5481DFA3B5E394F0.

![Image text](https://github.com/HwenZhang/TidalHydroFrac/blob/147148f5916b7197c94a07abe23951a49d448c2f/grounding_line_mesh_sensitivity/image/schematic.png)

# Governing equations

We solve for combined $\left(\boldsymbol u, p, \boldsymbol \tau\right)$, where $\boldsymbol{u}$ is the velocity, $p$ is the pressure, and $\boldsymbol\tau$ is the deviatoric stress. The function space is shown below.
```
P1 = FiniteElement('CG',triangle,1)     # Pressure
P2 = VectorElement('CG',triangle,2)     # Velocity
P3 = TensorElement("CG",triangle,2)    # Stress
```
## Mass and momentum conservation
Mass and momentum conservation gives
```math
\begin{gathered}
\boldsymbol{\nabla} \cdot \boldsymbol{u}=0 \\
\boldsymbol{\nabla} \cdot[-p \boldsymbol{I}+2 \eta(\boldsymbol{u}) \dot{\boldsymbol\varepsilon}]+\rho_i \boldsymbol{g}=0
\end{gathered}
```
where $\dot{\boldsymbol\varepsilon}$ is the strain rate, $\rho_i$ is the ice density and $\boldsymbol g$ is gravitational acceleration.

## Ice rheology: incompressible viscoelasticity
Ice is assumed to be incompressible upper-convective Maxwell-viscoelastic fluid. The constitutive law is
$$\boldsymbol{\tau}+\lambda\left(\frac{\partial \boldsymbol{\tau}}{\partial t}+\boldsymbol{u} \cdot \nabla \boldsymbol{\tau}-(\nabla \boldsymbol{u})^T \cdot \boldsymbol{\tau}-\boldsymbol{\tau} \cdot \nabla \boldsymbol{u}\right)-2 \eta(\boldsymbol{u}) \dot{\boldsymbol\varepsilon}=0$$, where $\lambda=\dfrac{\eta}{\mu}$ is the Maxwell time. In the computation, the characteristic Maxwell time is about $8~\text{hr}$.

### Nonlinear viscosity 
```math
\eta(\boldsymbol{u})=\frac{1}{2} B\left(|\boldsymbol{D}(\boldsymbol{u})|^2+\delta\right)^{-\frac{1}{3}},
```
where $B=2^{(n-1)/2n} A^{-1/n}$. In the computation, we use $A=3\times 10^{-24}~\text{Pa}^{-3}\cdot\text{s}^{-1}$, $n=3$. Note that $\delta=1\times 10^{-15}$ is the regularisation parameter to prevent infinite viscosity. 
### Elasticity
Shear modulus $\mu=0.18\times 10^{9}$ $\text{Pa}$.

# Boundary condition

## Left boundary 
Constant inflow
$$\boldsymbol u\cdot \hat{\boldsymbol x} = U_0 = 9\text{m/yr}$$ in x-direction.

## Right boundary
Ice overburden stress
$$\sigma_{xx}=-\rho_i g \left(h\left(x\right)-z\right),$$
representing purely floating ice shelf without internal deformation.

## Top and bottom boundary
Note in the computation, we use ALE method to catch the free surface movement. Free surface is governed by the kinematic equations
```math
\begin{gathered}
\frac{\partial h}{\partial t}(x, t)=\sqrt{1+\left(\frac{\partial h}{\partial x}\right)^2} u_n(x, h, t) \\
\frac{\partial s}{\partial t}(x, t)=-\sqrt{1+\left(\frac{\partial s}{\partial x}\right)^2} u_n(x, s, t),
\end{gathered}
```
where $h\left(x\right)$ and $s\left(x\right)$ represent the top and bottom surfaces of the ice shelf. The internal mesh is smoothed at each time step using
```
mesh.smooth()
```
.

## Bottom boundary
### Ice-bedrock interface
Penalty method to guarantee impenetrability. Friction is applied using a nonlinear sliding law
```math
\sigma_{nt}=-C\left[(\boldsymbol u\cdot \boldsymbol t)^2+\delta_s\right]^{-\frac{1}{3}}(\boldsymbol u\cdot \boldsymbol t),
```
where $C$ is a sliding constant which is set to measure the observed surface velocity. $\boldsymbol t$ is the tangential vector, $\delta_s=1\times 10^{-19}$ is the regularisation parameter. 

### Ice-ocean interface
Hydrostatic pressure, which is superposed by a base state sea level $l_0$ and a tidal perturbation, represented by a sinusoidal function.
```math
\sigma_{nn} = -\rho_w g \left[l_0 + A sin\left(2\pi f t\right)\right].
```

# Reference
Stubblefield, A. G.; Spiegelman, M.; Creyts, T. T. Variational Formulation of Marine Ice-Sheet and Subglacial-Lake Grounding-Line Dynamics. Journal of Fluid Mechanics 2021, 919, A23.
