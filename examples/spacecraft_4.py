"""Demo."""
import sympy as sp
import numpy as np
from f4cs.RWS_disturbances import RWSDisturbed
from f4cs.synthesis import TBS
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def saturate(x):
    """Shorthand for saturation function."""
    return sp.Max(-bound, sp.Min(bound, x))


# Define state and input variables
var_list = x1, x2, x3 = sp.symbols('x1,x2,x3')
input_list = u1, u2, u3 = sp.symbols('u1, u2, u3')
dist_var = w1, w2, w3 = sp.symbols('w1, w2, w3')
n = len(var_list)
m = len(input_list)

# Matrix/vector form.
xvec = sp.Matrix(var_list)
wvec = sp.Matrix(dist_var)
yvec = xvec+wvec
uvec = sp.Matrix(input_list)

# Inertia matrix and its inverse
inertia = inertia = sp.diag(4192., 342., 4200.)
inv_inertia = inertia.inv()

# Input saturation
bound = 50
u_sat = sp.Matrix([saturate(u1), saturate(u2), saturate(u3)])

# Define the system dynamics
f_sym = inv_inertia*(xvec.cross(inertia*xvec)) + \
    inv_inertia*u_sat

# Define the system specification sets and disturbance set
Slist = [[-2, 2], [-2, 2], [-2, 2]]
Ilist = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
Olist = [[-0.05, 0.05], [-0.05, 0.05], [-0.05, 0.05]]
Omega_list = [[-0.01, 0.01], [-0.01, 0.01], [-0.01, 0.01]]

# Assuming a quadratic LF and linear controller,
p = sp.symbols('p0:{}'.format(n**2+n*m+1))
p_k = p[0:n*m]
p_v = p[n*m:n**2+n*m]
x_vec = sp.Matrix(var_list)
K = sp.Matrix(sp.Array(p_k).reshape(m, n))
Q = sp.Matrix(sp.Array(p_v).reshape(n, n))

kappa = sp.flatten(sp.Array((K*yvec).T))
V = (x_vec.T*Q*x_vec)[0] + p[n**2+n*m]

# Define a controller and certificate function template
k_template = kappa
v_template = V

# Create a template dictionary. If 'constants' is kept empty, default values
# are used.
template = {'controller': k_template,
            'certificate': v_template,
            'parameters': p
            }

# path where the SMT files will be stored
path = 'e:/docker_connect/data'

# Use dReal and set its specific options
smt_options = {'solver': 'dReal', 'path': path, 'dprecision': 0.01}

options = {'Slist': Slist,  # Interval list of the safe set
           'Ilist': Ilist,  # Interval of the initial set
           'Olist': Olist,  # Interval list of the goal set
           'auxlist': Omega_list,
           'aux_var': dist_var,
           'number_samples': 50000,  # Number of (initial) samples
           # Maximum number of counterexamples
           'max_samp': 52000,
           'gamma': 0.01,  # (arbitrary) decrease of the LF
           'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
           'smt_options': smt_options,
           'epsilon': 0.1,  # Robustness buffer for the sample-based fitness
           'max_iterations': 30
           }

# Initiate the specification
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
