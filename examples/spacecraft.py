"""
Simple spacecraft example.

In this example, the angular velocities are stabilized, while never exceeding
a maximum angular velocity.

The system is subjected to input constraints and measurement disturbances.

To use, please revise 'path' and 'dReal_path' for your specific system.

@author: Cees F. Verdier
"""
import sympy as sp
import numpy as np
from f4cs.specifications.reach_while_stay import RWS
from f4cs.synthesis.template_based_synthesis import (
    TemplateBasedSynthesis as TBS)
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# CHANGE TO YOUR SYSTEM PATHS:

# path where the SMT files will be stored
path = 'e:/docker_connect/data'
# Path of dReal (Mac/Linux only)
dReal_path = '/opt/dreal/4.18.02.2/bin/dreal'


def saturate(x):
    """Shorthand for saturation function."""
    return sp.Max(-bound, sp.Min(bound, x))


# Define state and input variables
var_list = x1, x2, x3 = sp.symbols('x1,x2,x3')
input_list = u1, u2, u3 = sp.symbols('u1, u2, u3')
dist_var = w1, w2, w3 = sp.symbols('w1, w2, w3')

# Matrix/vector form.
x = sp.Matrix(var_list)
w = sp.Matrix(dist_var)
u = sp.Matrix(input_list)
y = x+w

# Inertia matrix and its inverse
J = sp.diag(4192., 342., 4200.)
J_inv = J.inv()

# Input saturation
bound = 50
u_sat = sp.Matrix([saturate(u1), saturate(u2), saturate(u3)])

# Define the system dynamics
f_sym = J_inv*(x.cross(J*x)) + J_inv*u_sat

# Define the system specification sets and disturbance set
S_list = [[-2, 2], [-2, 2], [-2, 2]]
I_list = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
O_list = [[-0.05, 0.05], [-0.05, 0.05], [-0.05, 0.05]]
Omega_list = [[-0.01, 0.01], [-0.01, 0.01], [-0.01, 0.01]]

# Assuming a quadratic LF and linear controller,
n = len(var_list)
m = len(input_list)
p = sp.symbols('p0:{}'.format(n**2+n*m+1))
K = sp.Matrix(p[0:n*m]).reshape(m, n)
Q = sp.Matrix(p[n*m:n**2+n*m]).reshape(n, n)

k_template = K*y
v_template = (x.T*Q*x)[0] + p[n**2+n*m]

# Create a template dictionary.
template = {'controller': k_template,
            'certificate': v_template,
            'parameters': p
            }

# Use dReal
smt_options = {'solver': 'dReal', 'path': path, 'dReal_precision': 0.001,
               'dReal_path': dReal_path,  # Path of dReal. Only for Mac/Linux
               't_max': 20  # Maximum time the SMT solver is running
               }

options = {'variables': var_list,  # List of symbolic states
           'inputs': input_list,  # List of symbolic inputs
           'uncertainties': dist_var,  # List of disturbances
           'system': f_sym,  # System dynamics
           'S_list': S_list,  # Interval list of the safe set
           'I_list': I_list,  # Interval of the initial set
           'O_list': O_list,  # Interval list of the goal set
           'uncertainty_list': Omega_list,  # Interval of the disturbances
           'number_samples': 50000,  # Number of (initial) samples
           'max_samp': 52000,  # Maximum number of counterexamples
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,  # SMT solver options
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 20,
           }

# Initiate the specification
spec = RWS(options)
# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)

# Repeated test
number_runs = 10
for i in range(number_runs):
    # Synthesize a solution.
    solution = synthesiser.synthesis(spec, template)

# Extract the controller (plus some input disturbances)
t = sp.symbols('t')
k = solution.k_sym.subs(zip(solution.p, solution.par))

# Run multiple simulations
tvec = np.linspace(0, 100, 2000)
number_trajectories = 10
y = [tvec]*number_trajectories

for i in range(number_trajectories):
    # random initial condition
    y0 = 0.5 * np.random.randint(-1, 2, size=3)
    # Random disturbances
    d_par = np.random.rand(7)
    disturbances = sp.Matrix([0.01*sp.sin(t*d_par[0]+d_par[1]),
                              0.01*sp.cos(t*d_par[2]+d_par[3]),
                              0.01*sp.sin(t*d_par[4]+d_par[5])])
    k = k.subs(zip(dist_var, disturbances))

    # Closed-loop system
    f = f_sym.subs(zip(input_list, k))
    f_fun = sp.lambdify([var_list, t], sp.flatten(f))

    y[i] = odeint(f_fun, y0, tvec)

# Plotting
y1 = np.array([x[:, 0] for x in y])
y2 = np.array([x[:, 1] for x in y])
y3 = np.array([x[:, 2] for x in y])

fig0, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
ax0.plot(tvec, y1.T, tvec, tvec*0+0.05, 'k--', tvec, tvec*0-0.05, 'k--')
ax0.set_ylabel('x1')
# ax0.set_ylim([-0.8, 0.8])
ax1.plot(tvec, y2.T, tvec, tvec*0+0.05, 'k--', tvec, tvec*0-0.05, 'k--')
ax1.set_ylabel('x2')
# ax1.set_ylim([-0.8, 0.8])
ax2.plot(tvec, y3.T, tvec, tvec*0+0.05, 'k--', tvec, tvec*0-0.05, 'k--')
ax2.set_ylabel('x3')
# ax2.set_ylim([-0.8, 0.8])
ax2.set_xlabel('time [s]')
plt.show()
