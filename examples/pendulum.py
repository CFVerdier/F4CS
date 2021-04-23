# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 13:58:59 2021

Simple pendulum example.

To use, please revise 'path' and 'dReal_path' for your specific system.

@author: Cees F. Verdier
"""
import sympy as sp
from f4cs.specifications.reach_while_stay import RWS
from f4cs.synthesis import TBS

# CHANGE TO YOUR SYSTEM PATHS:

# path where the SMT files will be stored
path = 'e:/docker_connect/data'
# Path of dReal (Mac/Linux only)
dReal_path = '/opt/dreal/4.18.02.2/bin/dreal'

# Define state and input variables
var_list = x1, x2 = sp.symbols('x1,x2')
input_list = u1, = sp.symbols('u1,')

# Define your Dynamics
f_sym = sp.Matrix([x2, 19.6*sp.sin(x1)-16*x2+4*sp.cos(x1)*u1])

Slist = [[-6, 6], [-10, 10]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5]]
Olist = [[-0.25, 0.25], [-0.25, 0.25]]

# Create a vector of constants for the template, as well as a symbolic variant
constants = [1, 1, 1, 1, 1, 1]
p = sp.symbols('p0:{}'.format(len(constants)))

# Define a controller and certificate function template
k_template = sp.Array([p[0]*x1+p[1]*x2])
v_template = p[2]*x1**2+p[3]*x2**2+p[4]*x1*x2+p[5]

# Create a template dictionary. If 'constants' is kept empty, default values
# are used.
template = {'controller': k_template,
            'certificate': v_template,
            'parameters': p,
            'values': constants  # optional. Otherwise, initialized with zeros
            }

# Use dReal
smt_options = {'solver': 'dReal', 'path': path, 'dprecision': 0.01,
               'dReal_path': dReal_path  # Path of dReal. Only for Mac/Linux
               }

options = {'Slist': Slist,  # Interval list of the safe set
           'Ilist': Ilist,  # Interval list of the initial set
           'Olist': Olist,  # Interval list of the goal set
           'number_samples': 10000,  # Number of (initial) samples
           # Maximum number of samples (when adding violations)
           'max_samp': 13000,
           'rdelta': 0.01,  # Inflation of the boundary
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }
# Initialize specification
spec = RWS(var_list, input_list, f_sym, options)
# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.

solution = synthesiser.synthesis(spec, template)
