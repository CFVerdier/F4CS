# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 13:58:59 2021

@author: ceesv
"""
import sympy as sp
from f4cs.specifications import RWS
from f4cs.synthesis import TBS

# Define state and input variables
var_list = x1, x2 = sp.symbols('x1,x2')
input_list = u1, = sp.symbols('u1,')

# Define your Dynamics
f_sym = sp.Matrix([x2, 19.6*sp.sin(x1)-16*x2+4*sp.cos(x1)*u1])

Slist = [[-6, 6], [-10, 10]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5]]
Olist = [[-0.25, 0.25], [-0.25, 0.25]]


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
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }
# Initialize specification
spec = RWS(var_list, input_list, f_sym, options)

synthesiser = TBS(options)
solution = synthesiser.synthesis(spec)


