# Solving-Non-linear-PDE
## Objective of the work
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
## Initial and Boundary Conditions
Initial condition is again in general form, and boundary condition two options are considered: Drichilet and Neumann, and it can be changed as required.
## Steps Considered
## Step-1
Write a finite difference based code and which used CN method to generate the data. To get the data following variable values are considered.

$$
v(x,t)=0.5\left(1+0.2\sin\left(\frac{2\pi t}{T}\right)\right)
$$

$$
D(x,t)=0.05\left(1+0.5\frac{x}{L}\right)
$$

$$
S(x,t)=0.02\exp\left(-\frac{(x-3)^2}{0.5}\right)
\begin{cases}
1, & t<5,\\
0, & t\geq 5.
\end{cases}
$$

