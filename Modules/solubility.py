# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 16:38:18 2021

@author: marcell
"""

from numpy import exp, log
from pyswarm import pso
from scipy.optimize import minimize


def func_solubility(T, Theta):
    return exp(Theta[0] + Theta[1] * T + Theta[2] * log(T))

def objective_function(Theta, T, b):
    f = func_solubility(T, Theta)
    sol = 0.0
    TF = len(f)
    for i in range(0,TF):
        sol = sol + (b[i] - f[i])**2
    return sol

def n_linear_regression(T, b):
    lb = [-20.0, -20.0, -20.0]
    ub = [20.0, 20.0, 20.0]
    
    #Heuristic Optimizer
    xopt, fopt = pso(objective_function, lb, ub, ieqcons=[], f_ieqcons=None, args=(T, b), kwargs={}, swarmsize=100,
                 omega=0.5, phip=0.8, phig=0.8, maxiter=5000, minstep=1e-6, minfunc=1e-6, debug=False)
    paramr = xopt
    x0 = (paramr)
    
    #Deterministic Optimizer
    Thetas = minimize(objective_function,x0,args=(T, b),method='Nelder-Mead', tol = 1e-7, options = {'xatol':1e-7, 'disp':True, 'maxiter': 2000})

    return Thetas.x
