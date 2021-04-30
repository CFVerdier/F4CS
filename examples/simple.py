"""

Simple  example.

@author: Cees F. Verdier
"""
import sympy as sp
from f4cs.specifications.reach_while_stay import RWS
from f4cs.synthesis.template_based_synthesis import (
    TemplateBasedSynthesis as TBS)
from f4cs.plotting import plotting

# Define state and input variables
var_list = x1, x2 = sp.symbols('x1,x2')
input_list = u1, = sp.symbols('u1,')

# Define your Dynamics
f_sym = sp.Matrix([x2, u1])  # Column vector

# RWS specification
S_list = [[-15, 15], [-15, 15]]
I_list = [[-5, 5], [-5, 5]]
O_list = [[-0.5, 0.5], [-0.5, 0.5]]

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

options = {'variables': var_list,  # tuple of symbolic variables
           'inputs': input_list,  # tuple of symbolic inputs
           'system': f_sym,  # symbolic equations of motion
           'S_list': S_list,  # Interval list of the safe set
           'I_list': I_list,  # Interval list of the initial set
           'O_list': O_list,  # Interval list of the goal set
           'number_samples': 1000,  # Number of (initial) samples
           # Maximum number of samples (when adding violations)
           'max_samples': 1300,
           'r_delta': 0.01,  # Inflation of the boundary
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }
# Initialize specification
spec = RWS(options)
# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.
# solution = synthesiser.synthesis(spec, template)

# Repeated test
number_runs = 1
for i in range(number_runs):
    # Synthesize a solution.
    solution = synthesiser.synthesis(spec, template)

# Extract controller
k = solution.k_sym.subs(zip(solution.p, solution.par))
f_closed_loop = f_sym.subs(zip(input_list, k))

# Plot the result
plotting(var_list, f_closed_loop, I_list)
