# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:26:48 2021

@author: ceesv
"""
import sympy as sp
from specifications import RWS
from candidates import Solution
import cma

# import timeit

var_list = x1, x2 = sp.symbols('x1,x2')
input_list = u1, = sp.symbols('u1,')

# Dynamics
f_sym = sp.Matrix([x2, u1])  # Column vector

Slist = [[-15, 15], [-15, 15]]
Ilist = [[-5, 5], [-5, 5]]
Olist = [[-1, 1], [-1, 1]]

# path where the SMT files will be stored
path = 'e:/docker_connect/data'

# Use Z3
# smt_options = {'solver': 'Z3'}

# Use dReal
smt_options = {'solver': 'dReal', 'path': path, 'dprecision': 0.01}

options = {'Slist': Slist,  # Interval list of the safe set
           'Ilist': Ilist,  # Interval list of the initial set
           'Olist': Olist,  # Interval list of the goal set
           'number_samples': 100,  # Number of (initial) samples
           # Maximum number of samples (when adding violations)
           'max_samp': 300,
           'rdelta': 0.01,  # Inflation of the boundary
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1}  # Robustness buffer for the sample-based fitness

# Initialize specification
spec = RWS(var_list, input_list, f_sym, options)
# Create an (hardcoded) individual
ind = Solution(spec)

# Give the individual arbitrary new parameters for testing
ind.par = [1, 1, 1]
ind.substitute_parameters(ind.par)

sigma0 = 0.5
cma.fmin(spec.parameter_fitness, ind.par, sigma0,
         args={ind, }, options={'verbose': -9})
spec.verify(ind)
print(spec.verification_result)
