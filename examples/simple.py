"""

Simple  example.

@author: Cees F. Verdier
"""
import sympy as sp
from f4cs.specifications import RWS
from f4cs.synthesis import TBS

# Define state and input variables
var_list = x1, x2 = sp.symbols('x1,x2')
input_list = u1, = sp.symbols('u1,')

# Define your Dynamics
f_sym = sp.Matrix([x2, u1])  # Column vector

# RWS specification
Slist = [[-15, 15], [-15, 15]]
Ilist = [[-5, 5], [-5, 5]]
Olist = [[-1, 1], [-1, 1]]

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

# Use Z3
smt_options = {'solver': 'Z3'}

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
# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.
solution = synthesiser.synthesis(spec, template)


