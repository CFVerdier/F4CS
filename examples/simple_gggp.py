"""

Simple  example.

@author: Cees F. Verdier
"""
import sympy as sp
from f4cs.specifications.reach_while_stay import RWS
from f4cs.synthesis.gggp.grammar_guided_genetic_programming import (
    GrammarGuidedGeneticProgramming as GGGP)
from f4cs.synthesis.gggp.grammar import Grammar
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

# Define the grammar
P = {'pol': ('pol+pol', 'const*mon'),
     'mon': ('mon*mon', 'vars'),
     'vars': ('x1', 'x2')
     }
P_terminal = {'pol': ('const*mon',),
              'mon': ('vars',),
              'vars': ('x1', 'x2')
              }
# Special production rules for constants
P_constants = {'const': ('Real(-1,1)',)}
# start = ('pol+const', sp.Matrix(['pol' for _ in range(0, len(input_list))]))
start = ('const*x1**2 + const*x2**2 + const*x1*x1 +const',
         sp.Matrix(['pol' for _ in range(0, len(input_list))]))
grammar = Grammar(start, P, P_terminal, P_constants)
max_depth = 5
number_individuals = 6

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
           'max_iterations': 30,
           'grammar': grammar,
           'max_depth': max_depth,
           'number_individuals': number_individuals
           }
# Initialize specification
spec = RWS(options)
# Initialize an template-based synthesis procedure.
synthesiser = GGGP(options)
print(synthesiser.population[0])
# Synthesize a solution.
solution = synthesiser.synthesis(spec)
print('optimized')
print(synthesiser.population[0])
print('fitness')
print(synthesiser.sample_fitness)

# # Repeated test
# number_runs = 1
# for i in range(number_runs):
#     # Synthesize a solution.
#     solution = synthesiser.synthesis(spec, template)

# # Extract controller
# k = solution.k_sym.subs(zip(solution.p, solution.par))
# f_closed_loop = f_sym.subs(zip(input_list, k))

# # Plot the result
# plotting(var_list, f_closed_loop, I_list)
