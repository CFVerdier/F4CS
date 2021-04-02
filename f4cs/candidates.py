# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:24:34 2021

@author: ceesv
"""

import sympy as sp
import numpy as np
from specifications import RWS

# TODO use the class below for new individuals
# TODO: use @property(?) for parameters?
# TODO: Remove spec specific data from the solution.


def custom_amax(x, **kwargs):
    """Alternative maximum for lambdafication."""
    return np.maximum(x[0], x[1])


def custom_amin(x, **kwargs):
    """Alternative minimum for lambdafication."""
    return np.minimum(x[0], x[1])


class Solution:
    """ Stores and manipulates a candidate solution for RWS.

    Constructor arguments:
        S (spec): specification from the RWS class

    Optimal arguments:
        none

    Attributes:
        TODO

    Methods:
        none
    """
    # TODO: variable in the entities that it can has: LF, controller, auxillary functions

    def __init__(self, S):

        if isinstance(S, RWS):  # Check if it is a RWS spec.
            # TODO replace with e.g. a dictionary to switch modes
            # HARDCODED CONTROLLER AND CLBF
            # self.par = [-1, -2, -78]

            self.spec = S

            self.par = [1, 1, 1]
            self.par_len = len(self.par)
            self.p = sp.symbols('p0:{}'.format(self.par_len))
            p = self.p
            x = self.spec.var
            self.k_sym_p = sp.Array([p[0]*x[0]+p[1]*x[1]])
            self.V_sym_p = x[0]**2+x[1]**2+x[0]*x[1]+p[2]

            self.k_sym = self.k_sym_p
            self.V_sym = self.V_sym_p
            self.dV_sym = sp.diff(self.V_sym, [self.spec.var])
            self.dtV_sym = sp.Matrix(self.dV_sym).dot(
                self.spec.f_sym).subs(zip(self.spec.input, self.k_sym))
            self.dV_sym = None
            self.dVfun = None

            # self.k_sym = None
            # self.V_sym = None
            # self.dV_sym = None
            # self.dVfun = None
            # self.k_fun = None
            # self.V_fun = None
            # self.dtV_sym = None

            self.fitness = None

            # Create the conditions dependent on the specification
            self.conditions_sym = self.spec.create_conditions(self)
            # Create symbolic fitness function for this individual.
            self.create_fitness()

            # Create conditions for the current parameters
            self.conditions = None
            self.substitute_parameters()

            # substitute the parameter vector and create the functions
            # self.substitute_parameters(self.par)
    def substitute_parameters(self):
        """Update conditions with current parameters."""
        self.conditions = tuple([condition.subs(
            zip(self.p, self.par)) for condition in self.conditions_sym])

    def create_fitness(self):
        """Create a symbolic fitness function given a symbolic condition."""
        result = ()
        for condition in self.conditions_sym:
            # Copy condition
            f = condition
            # Replace Or and And to prevent TypeError during the transformation
            f = f.replace(sp.Or, sp.Function('dummy_or'))
            f = f.replace(sp.And, sp.Function('dummy_and'))
            # Write to standard form.
            f = f.replace(sp.StrictGreaterThan, lambda x,
                          y: x >= y + self.spec.c)
            f = f.replace(sp.StrictLessThan, lambda x, y: x <= y - self.spec.c)
            f = f.replace(sp.LessThan, lambda x, y: -x >= -y)
            f = f.replace(sp.GreaterThan, lambda x, y: x-y-self.spec.epsilon)
            # Replace booleans
            f = f.replace(True, 0.0)
            f = f.replace(False, -np.inf)
            # Replace Or and And
            f = f.replace(sp.Function('dummy_or'), sp.Max)
            f = f.replace(sp.Function('dummy_and'), sp.Min)
            # Saturate minimal error
            f = sp.Min(f, 0.0)
            # Concatinate to result tuple
            function = sp.lambdify(
                self.p + self.spec.var, f,
                modules=[{'amax': custom_amax, 'amin': custom_amin}, "numpy"])
            result = result + (function,)
        self.fitness = result

    # def substitute_parameters(self, z):
    #     self.par = z
    #     self.k_sym = self.k_sym_p.subs(
    #         [(self.p[i], self.par[i]) for i in range(self.par_len)])
    #     self.V_sym = self.V_sym_p.subs(
    #         [(self.p[i], self.par[i]) for i in range(self.par_len)])
    #     self.dV_sym = sp.diff(self.V_sym, [self.S.var])
    #     self.dtV_sym = sp.Matrix(self.dV_sym).dot(
    #         self.S.f_sym).subs(zip(self.S.input, self.k_sym))

        # self.make_functions()

    # def make_functions(self):
    #     # Make a Controller and LBF function
    #     self.k_fun = sp.lambdify([self.S.var], self.k_sym, "numpy")
    #     self.V_fun = sp.lambdify([self.S.var], self.V_sym, "numpy")

    #     # Compute the derivatives and make into functions
    #     self.dV_sym = sp.diff(self.V_sym, [self.S.var])
    #     self.dV_fun = sp.lambdify([self.S.var], self.dV_sym, "numpy")
    #     self.dtV_sym = sp.Matrix(self.dV_sym).dot(
    #         self.S.f_sym).subs(zip(self.S.input, self.k_sym))
