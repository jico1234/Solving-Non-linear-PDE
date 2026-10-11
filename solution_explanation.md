## Solution Explanation
Non-linear ADE has been considered and solved numerically using CN method, which helps us to generate the data; after generation of data, PySR 
has been used to identify the relation between the data, initial condition, dispersion, velocity, and all the other functions are also kept with
the solution. Now, as we know that any solution must satisfy the equation, we input every solution to the equation and check the residue. Once we obtained
the residue is near zero and the R^2 value is more than .99, then we are saying that it is a solution of the PDE.
