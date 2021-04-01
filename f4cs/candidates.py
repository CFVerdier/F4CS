# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:24:34 2021

@author: ceesv
"""

import sympy as sp
from specifications import RWS

# TODO use the class below for new individuals


class Solution:
    """ Stores and manipulates a candidate solution for RWS

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
            self.S = S
            self.par = [-1, -2, -78]
            self.par_len = len(self.par)
            self.p = sp.symbols('p0:{}'.format(self.par_len))
            p = self.p
            x = S.var
            self.k_sym_p = sp.Array([p[0]*x[0]+p[1]*x[1]])
            self.V_sym_p = x[0]**2+x[1]**2+x[0]*x[1]+p[2]

            self.k_sym = None
            self.V_sym = None
            self.k_fun = None
            self.V_fun = None
            self.dV_sym = None
            self.dVfun = None
            self.dtV_sym = None

            # TODO: store these values automatically
            self.sample_fitness = 0
            self.fitness = 0

            # substitute the parameter vector and create the functions
            self.substitute_parameters(self.par)

            # TODO: use @property(?) for parameters?
            # TODO: Remove spec specific data from the solution.

    def substitute_parameters(self, z):
        self.par = z
        self.k_sym = self.k_sym_p.subs(
            [(self.p[i], self.par[i]) for i in range(self.par_len)])
        self.V_sym = self.V_sym_p.subs(
            [(self.p[i], self.par[i]) for i in range(self.par_len)])
        self.make_functions()

    def make_functions(self):
        # Make a Controller and LBF function
        self.k_fun = sp.lambdify([self.S.var], self.k_sym, "numpy")
        self.V_fun = sp.lambdify([self.S.var], self.V_sym, "numpy")

        # Compute the derivatives and make into functions
        self.dV_sym = sp.diff(self.V_sym, [self.S.var])
        self.dV_fun = sp.lambdify([self.S.var], self.dV_sym, "numpy")
        self.dtV_sym = sp.Matrix(self.dV_sym).dot(
            self.S.f_sym).subs(zip(self.S.input, self.k_sym))
