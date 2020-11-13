# Robust Quadratic Optimization 

This project contains all the data and codes related to the paper:
Marandi, Ahmadreza, et al. "Extending the scope of robust quadratic optimization." arXiv preprint arXiv:1909.01762 (2019).

## Table of contents
* [Technologies](#technologies)
* [Robust Portfolio Optimization] (#portfolio_optimization)
* [Robust Norm Approximation] (#norm_approximation)
* [Robust Regression Line] (#regression_line)

## Technologies
For this project, we use
* Matlab (the codes are written in Matlab R2019b)
* Yalmip (https://yalmip.github.io/) to pass optimization to solvers
* Mosek (https://www.mosek.com/) to solve the optimization problems

## Robust Portfolio Optimization
To find the mathematical formulation of the problem, see Section 7.1 of the paper.
To find a robust protfolio, we first need to construct an uncertainty set based on the historical data. We use the monthly average value weighted return of 5 and 30 industries from 1956 until 2015, obtained from  ``Industry Portfolios" data on the website http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html . For 5 industries, the data is stored in Data_5.mat . For 30 industries, the data is stored in Data_30.mat . The function ```$ RobustPortfolio.m ``` gets three inputs:
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


## Robust Regression Line
