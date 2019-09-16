# 2Level-HD
Implementation of a global search method for the two-level $\ell_1$ function in Matlab.
The two-level $\ell_1$ function is a nonconvex but piecewise linear sparse penalty.
Improving the hill detouring method, this global search method utilizes its piecewise linearity to escape from the local optima and can be used in convex feasible set.

# Main Files
- sim_optimization.m: a simulation experiment using GA, ABCD and this method to compare the optimization performance for the following problem.

       min_{x}     L_{two-level}(x) + (\|AX-b\|_2^2)/2*lambda        subject to l_i <= x_i <= u_i

- demo.m: a simulation experiment for noise-free compressive sensing problem. This method is compared with BP and ISD. You need YALL1 and ISD to call this file.

- solve_2LGD_noisefree.m: the descend method for solving the two-level $\ell_1$ function 

- solve_2LHD_noisefree.m: the global search method for solving the two-level $\ell_1$ function 

Both of these two function solving the following optimization problem.

       min_{x}     L_{two-level}(x)         subject to Phi*x=Y
