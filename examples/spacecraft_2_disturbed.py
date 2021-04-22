"""Demo."""
import sympy as sp
import numpy as np
from f4cs.RWS_disturbances import RWSDisturbed
from f4cs.synthesis import TBS
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define state and input variables
var_list = x1, x2, x3 = sp.symbols('x1,x2,x3')
input_list = u1, u2, u3 = sp.symbols('u1, u2, u3')
aux_var = w1, w2, w3 = sp.symbols('w1, w2, w3')


xvec = sp.Matrix(var_list)
uvec = sp.Matrix(input_list)
inertia = sp.Matrix([[4192., 0, 0], [0, 342., 0], [0, 0, 4200.]])
inv_inertia = inertia.inv()

bound = 50


def saturate(x):
    """Shorthand for saturation function."""
    return sp.Max(-bound, sp.Min(bound, x))


u_sat = sp.Matrix([saturate(u1), saturate(u2), saturate(u3)])

disturbances = sp.Matrix(aux_var)

# f_sym = inv_inertia*(xvec.cross(inertia*xvec)) + \
#     inv_inertia*u_sat + inv_inertia*disturbances

f_sym = inv_inertia*(xvec.cross(inertia*xvec)) + \
    inv_inertia*u_sat


Slist = [[-2, 2], [-2, 2], [-2, 2]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
Olist = [[-0.05, 0.05], [-0.05, 0.05], [-0.05, 0.05]]
auxlist = [[-0.01, 0.01], [-0.01, 0.01], [-0.01, 0.01]]


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
# kappa = sp.Array((K*x_vec).T)[0]

# Input disturbance
kappa = sp.flatten(sp.Array((K*(x_vec+disturbances)).T))

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
           'auxlist': auxlist,
           'aux_var': aux_var,
           'number_samples': 50000,  # Number of (initial) samples
           # Maximum number of counterexamples
           'max_samp': 52000,
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }

spec = RWSDisturbed(var_list, input_list, f_sym, options)


# Initialize an template-based synthesis procedure.
synthesiser = TBS(options)
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
    k = k.subs(zip(aux_var, disturbances))

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
