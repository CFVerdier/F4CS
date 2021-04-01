# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 21:35:41 2021

@author: ceesv
"""

import numpy as np
import sympy as sp
from specifications import Spec
import pyibex

#TODO: update to new specification format

class LocalStability(Spec):
    """Local stability specification.

    For general nonlinear systems, this class provides a method to abstract the
    system to a polynomial one, such that the conditions to verify become
    decidable.
    """

    def __init__(self, variables, inputs, f_sym, options):
        # Call the __init__ function of the Spec parent class first.
        Spec.__init__(self, variables, inputs, f_sym, options)

        self._number_conditions = 2  # number of local stability conditions

        D_list = self.options["Dlist"]
        # decrease of the LBF. Default 0.01
        self.gamma = self.options.get("gamma", 0.01)
        self.positive_definite = self.options.get(
            "PD", self.gamma*np.sum([xi**2 for xi in self.var]))

        # Create sample sets
        D_data = self.sample_set(D_list)
        self.data_sets = [D_data]

        # Create symbolic domains for SMT solver
        D_set = sp.And()
        for i in range(0, self.n):
            D_set = sp.And(
                D_set, sp.And(self.var[i] >= D_list[i][0],
                              self.var[i] <= D_list[i][1])
            )

        self.condition_set = (D_set, D_set)
        self.conditions = None
        self.verification_result = [None] * self._number_conditions

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        con1 = solution.V_sym >= self.positive_definite
        con2 = solution.dtV_sym <= -self.gamma*solution.V_sym

        self.conditions = (con1, con2,)

    def abstract_system(self, solution):
        """Abstract a general nonlinear system to a polynomial system."""
        f = solution.f_sym
        df = sp.diff(f, [self.var])
        ddf = sp.diff(df, [self.var])


# demo
var_list = x1, x2 = sp.symbols('x1,x2')
f = sp.Matrix([x2, sp.sin(x1)-0.1*x2+0.5*sp.cos(x1)*(-2*x1-0.1*x2)])
df = sp.diff(f, [var_list])
ddf = sp.diff(df, [var_list])

interval = pyibex.IntervalVector([[-1, 1], [-1, 1]])
my_fun = pyibex.Function('x1', 'x2', str(ddf[0, 0, 1, 0]))
my_fun.eval(interval)

n = len(var_list)

vec = [[0]*n for _ in range(n)]
remainder = [[[0]*n for _ in range(n)] for _ in range(n)]
symbolic_matrix = [[[0]*n for _ in range(n)] for _ in range(n)]
myfun = [0]*n
string_variables = [str(x) for x in var_list]
counter = 1
aux_list = []
aux_set_list = []
result_list = [[0]*n for _ in range(n)]
for i in range(0, n):
    for j in range(0, n):
        vec[i][j] = var_list[i]*var_list[j]
        ddfij = sp.diff(f, var_list[i], var_list[j])
        for k in range(0, n):
            # print(str(ddfij[k]))
            my_fun = pyibex.Function(*string_variables, str(ddfij[k]))
            a = my_fun.eval(interval)
            # print(a)
            remainder[i][j][k] = a
            if a.mag() == 0.0:
                symbolic_matrix[i][j][k] = 0
            else:
                new_aux = sp.Symbol('z'+str(counter))
                aux_set_list.append([a[0], a[1]])
                aux_list.append(new_aux)
                symbolic_matrix[i][j][k] = new_aux
                counter += 1
        result_list[i][j] = sp.Matrix(symbolic_matrix[i][j])*vec[i][j]

flat_results = [val for sublist in result_list for val in sublist]
result = sp.zeros(n, 1)
for element in flat_results:
    result = result+element
result = (0.5*result)
result.simplify
aux_list = tuple(aux_list)

point = sp.zeros(2,1)

# f_point = f.subs([(var_list[i], point[i]) for i in range(n)])
# df_point = df.subs([(var_list[i], point[i]) for i in range(n)])

# abstracted_f = f_sym(sp.zeros(2,1))

print(str(aux_list) + ' in ' + str(aux_set_list))
print(result)

# TODO: clean up code
# TODO: interface with Z3, clean up the verification class, etc
# TODO: actually synthesize a LF, hihi