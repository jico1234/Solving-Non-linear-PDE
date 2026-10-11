# Semi analytical soultion of Non-linear-PDE
## Objective of the work
Produce a closed solution which is semi-analytical in nature for nonlinear advection-diffusion equations. Where velocity, dispersion, source term, initial condition, and two boundary conditions can be considered as function of space and time. Can able to solve a family of PDE.
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
Write a finite difference based code and which used spatial discretization and an implicit time-integration scheme to generate the data. To get the data following variable values are considered.

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

$$
k=0.5
$$

Initial Condition

$$
C(x,0)=0
$$

$$
C(0,t)=1
$$

$$
\left.\frac{\partial C}{\partial x}\right|_{x=L}=0
$$

## Step-2
Used PySR to read the ade_data.csv file and remove al  t=0 and train on 3000 random rows and tests on 2000 random columns. It runs PySR with the operators +, -, *, / and the with 5 different unary function sets.Then finds the equations and store them to pysr_equations_1.csv to pysr_equations_5.csv consider all the values and check two propety residual<0.001 and R2 > .98. If both conditions hold, then we will consider that as a semi-analytical solution. If no such function has been found, then we will change our unary functions sets again and run the code.

## This is the Nobel procedure to find the semi-analytical solution using finite differences and PySR for nonlinear ADE.

## Limitations
1. With the increasing niterations, computation time may increase.
2. Some trivial boundary conditions are considered here; complex boundary conditions may take more computation time.
3. This is a one-dimensional problem; increase of dimension will use the same method but may need an additional linear solver.
