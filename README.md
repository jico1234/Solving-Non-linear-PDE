# Solving-Non-linear-PDE
**Objective of the work**
Produce a closed solution from a semi-analytical method for nonlinear advection-diffusion equations. It will feature variable speed, a source term for dispersion, the beginning condition, and two boundary conditions. Can able to solve a family of PDE.
## Governing Partial Differential Equation

The nonlinear advection–dispersion–reaction equation is given by

$$
\frac{\partial C}{\partial t} +\frac{\partial (vC)}{\partial x} = \frac{\partial}{\partial x} \left(D\frac{\partial C}{\partial x}\right) -kC^2+S
$$

where:

- $C(x,t)$: Contaminant concentration
- $v(x,t)$: Advection velocity
- $D(x,t)$: Dispersion coefficient
- $k$: Nonlinear reaction-rate coefficient
- $S(x,t)$: Source or sink term
- $x$: Spatial coordinate
- $t$: Time
