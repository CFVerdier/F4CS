"""Demo."""
import sympy as sp
import numpy as np
from f4cs.tracking import Tracking
from f4cs.specifications import RWS
from f4cs.synthesis import TBS

# Define state and input variables
var_list = x1, x2, x3, x4 = sp.symbols('x1,x2,x3,x4')
input_list = u1, u2, = sp.symbols('u1,u2')
aux_list = t, = sp.symbols('t,')

# Define your Dynamics: car example
# f = sp.Matrix([u1, u2, x1*sp.cos(x2), x1*sp.sin(x2)])
# x_ref = sp.Matrix([20+0.020567*t,
#                    0.19954*t,
#                    19.981 * t-0.10838*t**4,
#                    1.9915*t**2])
# u_ref = sp.Matrix([0.01835, 0.1995])

# xd_ref = sp.diff(x_ref, t)

# # Substitute into the old dynamics, to obtain the error dynamics
# f_sym = f.subs(zip(var_list, sp.Matrix(var_list)+x_ref))
# f_sym = f_sym.subs(zip(input_list, sp.Matrix(input_list)+u_ref))
# f_sym = f_sym - xd_ref

f_sym = sp.Matrix([x2, x3, x4, u1])

t_max = 10

# RWS specification
# Slist = [[-0.4, 0.4], [-0.2, 0.2], [-5, 5], [-5, 5]]
# Tlist = [[0, t_max]]
# Ilist = [[-0.15, 0.15], [-0.02, 0.02], [-0.2, 0.2], [-0.2, 0.2]]
# Olist = [[-0.1, 0.1], [-0.02, 0.02], [-0.2, 0.2], [-0.2, 0.2]]

# Slist = [[-100, 100], [-100, 100], [-100, 100], [-100, 100]]
# Ilist = [[-0.05, 0.05], [-0.01, 0.01], [-0.02, 0.02], [-0.02, 0.02]]
# Olist = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
# Tlist = [[0, t_max]]

Slist = [[-100, 100], [-100, 100], [-100, 100], [-100, 100]]
Tlist = [[0, t_max]]
Ilist = [[-0.15, 0.15], [-0.15, 0.15], [-0.15, 0.15], [-0.15, 0.15]]
Olist = [[-0.15, 0.15], [-0.15, 0.15], [-0.15, 0.15], [-0.15, 0.15]]


# Assuming a quadratic LF and linear controller
# Create a vector of constants for the template, as well as a symbolic variant
n = len(var_list)
m = len(input_list)

p = sp.symbols('p0:{}'.format(n**2+n*m+1))
constants = np.zeros(n**2+n*m+1)
p_k = p[0:n*m]
p_v = p[n*m:n**2+n*m]
x_vec = sp.Matrix(var_list)
K = sp.Matrix(sp.Array(p_k).reshape(m, n))
Q = sp.Matrix(sp.Array(p_v).reshape(n, n))
kappa = sp.Array((K*x_vec).T)[0]
V = (x_vec.T*Q*x_vec)[0] + p[n**2+n*m]

# Define a controller and certificate function template
k_template = kappa
v_template = V

# Create a template dictionary. If 'constants' is kept empty, default values
# are used.
template = {'controller': k_template,
            'certificate': v_template,
            'parameters': p,
            'values': constants  # optional. Otherwise, initialized with zeros
            }

# path where the SMT files will be stored
path = 'e:/docker_connect/data'

# Use dReal
smt_options = {'solver': 'dReal', 'path': path, 'dprecision': 0.01}

# Use Z3
# smt_options = {'solver': 'Z3'}

options = {'Slist': Slist,  # Interval list of the safe set
           'Ilist': Ilist,  # Interval list of the initial set
           'Olist': Olist,  # Interval list of the goal set
           'Tlist': Tlist,  # Interval list of the time interval
           'number_samples': 50000,  # Number of (initial) samples
           # Maximum number of samples (when adding violations)
           'max_samp': 52000,
           'rdelta': 0.01,  # Inflation of the boundary
           'alpha': 10,  # sublevel set constant
           'beta': 1,  # sublevel set constant for the goal set
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }
# Initialize specification
# spec = Tracking(var_list, input_list, f_sym, options, aux_list)
spec = RWS(var_list, input_list, f_sym, options)

# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.
solution = synthesiser.synthesis(spec, template)
