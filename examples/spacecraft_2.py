"""Demo."""
import sympy as sp
import numpy as np
from f4cs.specifications import RWS
from f4cs.practical_stability import PracticalStability
from f4cs.synthesis import TBS
import matplotlib.pyplot as plt

# Define state and input variables
var_list = x1, x2, x3 = sp.symbols('x1,x2,x3')
input_list = u1, u2, u3 = sp.symbols('u1, u2, u3')

xvec = sp.Matrix(var_list)
uvec = sp.Matrix(input_list)
inertia = sp.Matrix([[4192., 0, 0], [0, 342., 0], [0, 0, 4200.]])
inv_inertia = inertia.inv()

bound = 50


def saturate(x):
    """Shorthand for saturation function."""
    return sp.Max(-bound, sp.Min(bound, x))


u_sat = sp.Matrix([saturate(u1), saturate(u2), saturate(u3)])


f_sym = inv_inertia*(xvec.cross(inertia*xvec))+inv_inertia*u_sat

Slist = [[-2, 2], [-2, 2], [-2, 2]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
Olist = [[-0.05, 0.05], [-0.05, 0.05], [-0.05, 0.05]]


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
           'Ilist': Ilist,  # Interval of the initial set
           'Olist': Olist,  # Interval list of the goal set
           'number_samples': 50000,  # Number of (initial) samples
           # Maximum number of samples (when adding violations)
           'max_samp': 52000,
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }

spec = RWS(var_list, input_list, f_sym, options)


# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
# Synthesize a solution.
solution = synthesiser.synthesis(spec, template)

# fig = plt.figure(1)
# plt.scatter(spec.data_sets[1][:, 0], spec.data_sets[2][:, 1])