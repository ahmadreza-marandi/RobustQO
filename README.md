# Robust Quadratic Optimization 

This project contains all the data and codes related to the paper:
* Marandi, Ahmadreza, et al. "Extending the scope of robust quadratic optimization." arXiv preprint arXiv:1909.01762 (2019). 
In case of using these codes, you should properly cite the paper.  

## Table of contents
* [Technologies](#technologies)
* [Robust Portfolio Optimization](#Robust-Portfolio-Optimization)
* [Robust Norm Approximation](#robust-norm-approximation)
* [Robust Regression Line](#robust-regression-line)

## Technologies
For this project, we use
* Matlab (the codes are written in Matlab R2019b)
* Yalmip (https://yalmip.github.io/) to pass optimization to solvers
* Mosek (https://www.mosek.com/) to solve the optimization problems

## Robust Portfolio Optimization
To find the mathematical formulation of the problem, see Section 7.1 of the paper.
To find a robust protfolio, we first need to construct an uncertainty set based on the historical data. We use the monthly average value weighted return of 5 and 30 industries from 1956 until 2015, obtained from  "Industry Portfolios" data on the website http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html . For 5 industries, the data is stored in Data_5.mat . For 30 industries, the data is stored in Data_30.mat . The function ```$ RobustPortfolio.m ``` gets three inputs:
* z, which is the historical data
* alpha, which is linked to the size of the uncertainty set (see the discussion on Section 7.1 of the paper)
* lambda, which is the risk-aversion coefficient.
Then the optimization problem is constructed and gets solved. The function then returns five outputs:
* z_a, which is the probability that the objective function is higher than the obtained optimal value
* Rphi, which is the porfolio selection solution
* W, which is the optimal solution for matrix W added in the reformulation (see problem (30) in the paper)
* s_moment, which is the size of the moment matrix V (see (21) in the paper)
* solver_time, which is the time taken by the solver to find the optimal solution

We emphasize that for this function, function ``` $svec.m``` is needed. 

Next to obtaining the robust solution, you can also obtain a stochastic solution (solution of optimization problem containing chance constraint). Function ``` $ ChancePortfolio.m ``` is the function to so (see equation (26) in the paper). This function gets the same inputs as ``` $RobustPortfolio.m ``` but returns
* Chancephi, which is the portfolio solution,
* ChanceW, which is the part of optimal solution associated to the matrix W (see equation (27))
* s_moment, which is the size of the moment matrix V
* R, which is the Cholesky factorization of the moment matrix V
* solver_time, which is the time taken by the solver to find the solution.


An example of how we use these functions is provided in ``` $ Portfolio_30.m ``` .




## Robust Norm Approximation
For norm approximation problems, we work with ill-conditioned matrices. The main reason for that is that for ill-conditioned matrices, the nominal solution can be extremly sensitive to an error to the data. Hence, we generate 10 random instances with the size of n, varying between 4 until 100. Moreover, we have construct different uncertainty sets, which can be seen as a generalized budget uncertainty set. This data is stored in ``` Random_instances.mat ``` . In this file, you will find three cells: 

* matrix, with dimension 10x100. The first dimension is for the index of the random instance and the second dimension is for the size of the matrix A
* right, with dimension 10x100. The first dimension is for the index of the random instance and the second dimension is for the size of the matrix b
* Uncertainty, with dimension 10x100x3. The first dimension is for the index of the random instance, the second dimension is for the size, and the third one contains B1, B2, and K:

```
B1=Uncenrtainty{random_construction,n,1}
B2=Uncenrtainty{random_construction,n,2}
K=Uncenrtainty{random_construction,n,3};
A=matrix{random_construction,n};
b=right{random_construction,n};
 ``` 
For formulations, see Section 7.2. 

These inputs are used to check the quality of the method. The codes related to optimization is provided in ``` Norm_approximation_randomUncertain.m ``` . Using this file, you can solve the inner and outer approximations (equations (32) and (33)). To evaluate the solutions, we developed a heuristic method to find the worst-case scenarios for a given solution. The function ``` worst_value.m ``` coded this heuristic method. This function has eight inputs:

* rho, related to the size of the uncertainty set
* B1, the coefficient used in constructing the uncertainty set,
* B2, another coefficient used constructing the uncertainty set,
* y, is the solution we are analyzing, 
* x_value, is the value of Ay-b,
* K, is linked to the size of the uncertainty set,
* m, is the number of rows in A
* n, is the size of the solution.

For the description of the heuristic method, see Appendix D.

## Robust Regression Line

This folder contains the data and the code related to solving a robust linear regression problem. We use the data of the papar:

Candanedo, Luis M., Véronique Feldheim, and Dominique Deramaix. "Data driven prediction models of energy use of appliances in a low-energy house." Energy and buildings 140 (2017): 81-97.

The code simply solves the optimization problem in Appendix E. The ```worst_value.m ``` function is the same function, explained above, to find an approximation of the worst-case scenario given a solution. 
