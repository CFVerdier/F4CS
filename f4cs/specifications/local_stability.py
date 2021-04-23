# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 21:35:41 2021

EXPERIMENTAL MODULE

@author: Cees F. Verdier
"""

import numpy as np
import sympy as sp
from .specification import Specification
import pyibex

# TODO: update to new specification format

print("Warning: this is an experimental module")

class LocalStability(Specification):
    """Local stability specification.

    For general nonlinear systems, this class provides a method to abstract the
    system to a polynomial one, such that the conditions to verify become
    decidable.
    """

    def __init__(self, variables, inputs, f_sym, options):
        # Call the __init__ function of the Spec parent class first.
        Specification.__init__(self, variables, inputs, f_sym, options)

        self._number_conditions = 2  # number of local stability conditions

        D_list = self.options["Dlist"]
        self.x0 = self.options.get('x0', sp.zeros(self.n, 1))
        # decrease of the LBF. Default 0.01
        self.gamma = self.options.get("gamma", 0.01)
        self.positive_definite = self.options.get(
            "PD", self.gamma*np.sum([xi**2 for xi in self.var]))

        # Create sample sets
        D_data = self.sample_set(D_list)
        self.data_sets = [D_data, D_data]

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

        self.f_abstracted = None
        # Rewrite the system with the abstracted system
        self.abstract_system()

    def abstract_system(self):
        """Abstract a general nonlinear system to a polynomial system."""
        f = self.f_sym
        interval = pyibex.IntervalVector(self.options["Dlist"])
        x = sp.Matrix(self.var)
        f0 = f.subs(zip(self.var, self.x0))
        df = f.jacobian(x)
        f1 = (df.subs(zip(self.var, self.x0)))*(x-self.x0)
        symbolic_sub_matrix = [0]*self.n
        counter = 1
        string_variables = [str(x) for x in self.var]
        aux_list = ()
        aux_set_list = []
        f2 = sp.zeros(self.n, 1)
        for i in range(0, self.n):
            for j in range(0, self.n):
                vec = (self.var[i]-self.x0[i])*(self.var[j]-self.x0[j])
                ddfij = sp.diff(f, self.var[i], self.var[j])
                for k in range(0, self.n):
                    my_fun = pyibex.Function(*string_variables, str(ddfij[k]))
                    a = my_fun.eval(interval)
                    if a.mag() == 0.0:
                        # No interval at this position.
                        symbolic_sub_matrix[k] = 0
                    else:
                        # Make a new variable
                        new_aux = sp.Symbol('z'+str(counter))
                        aux_set_list.append([a[0], a[1]])
                        aux_list = aux_list + (new_aux,)
                        # Add this variable at this position in the matrix
                        symbolic_sub_matrix[k] = new_aux
                        counter += 1
                f2 = f2 + 0.5*sp.Matrix(symbolic_sub_matrix)*vec

        # create sets
        # TODO: automate this function
        aux_set = sp.And()
        for i in range(len(aux_set_list)):
            aux_set = sp.And(
                aux_set, sp.And(aux_list[i] >= aux_set_list[i][0],
                                aux_list[i] <= aux_set_list[i][1])
            )
        # TODO: hackish

        aux_data = self.sample_set(aux_set_list, self.number_samples)
        for i in range(self._number_conditions):
            self.data_sets[i] = np.append(self.data_sets[i],
                                          aux_data, axis=1)

        self.f_sym = f0 + f1 + f2
        self.var = self.var + aux_list
        new_conditions = [sp.And(condition, aux_set)
                          for condition in self.condition_set]
        self.condition_set = tuple(new_conditions)

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        # con1 = solution.V_sym >= self.gamma*self.positive_definite
        # con2 = solution.dtV_sym <= -self.gamma*solution.V_sym
        con1 = solution.V_sym >= self.positive_definite
        con2 = solution.dtV_sym <= -self.positive_definite

        return (con1, con2,)

    def demo():
        """Short demo of the class."""
        from candidates import Solution
        import cma

        var_list = x1, x2 = sp.symbols('x1,x2')

        f_sym = sp.Matrix([x2,
                           19.6*sp.sin(x1)-16*x2 +
                           4*sp.cos(x1)*(2.624*x1+1.0516*x2)])

        D_list = [[-0.25, 0.25], [-0.25, 0.25]]

        smt_options = {'solver': 'Z3'}


        options = {'Dlist': D_list,  # Interval list of the domain
                   'number_samples': 100,  # Number of (initial) samples
                   # Maximum number of samples (when adding violations)
                   'max_samp': 300,
                   'c': 0.01,  # (arbitrary) nonnegative parameter (see manual)
                   'smt_options': smt_options,
                   'epsilon': 0.1}  # Robustness buffer of sample-based fitness

        spec = LocalStability(var_list, (), f_sym, options)
        # Create an (hardcoded) individual
        ind = Solution(spec)

        sigma0 = 0.5
        cma.fmin(spec.parameter_fitness, ind.par, sigma0,
                 args={ind, }, options={'verbose': -1})

        spec.verify(ind)
        print(spec.verification_result)
