# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 13:58:59 2021

@author: ceesv
"""
import sympy as sp
import numpy as np
from f4cs.specifications import RWS
from f4cs.synthesis import TBS

# Define state and input variables
var_list = x1, x2, x3 = sp.symbols('x1, x2, x3')
input_list = u1, u2, = sp.symbols('u1, u2')

# Define your Dynamics
f_sym = sp.Matrix([u1, u2, x1*x2])

Slist = [[-5, 5], [-5, 5], [-5, 5]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5], [1, 2]]
Olist = [[-0.2, 0.2], [-0.2, 0.2], [-0.2, 0.2]]

# Create a vector of constants for the template, as well as a symbolic variant
constants = np.ones(13)
p = sp.symbols('p0:{}'.format(len(constants)))

# Define a controller and certificate function template
k_template = sp.Array([p[0]*x1+p[1]*x2+p[2]*x3, p[3]*x1+p[4]*x2+p[5]*x3])
v_template = p[6]*x1**2+p[7]*x2**2+p[8] * \
    x3**2+p[9]*x1*x2+p[10]*x1*x3+p[11]*x2*x3+p[12]

# Create a template dictionary. If 'constants' is kept empty, default values
# are used.
template = {'controller': k_template,
            'certificate': v_template,
            'parameters': p,
            'values': constants  # optional. Otherwise, initialized with zeros
            }


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
           'max_iterations': 100
           }
# Initialize specification
spec = RWS(var_list, input_list, f_sym, options)

# Initialize specification
spec = RWS(var_list, input_list, f_sym, options)

# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.
solution = synthesiser.synthesis(spec, template)
