# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:24:34 2021

@author: Cees F. Verdier
"""

import sympy as sp
import numpy as np

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
    """Stores and manipulates a candidate solution for RWS.

    Constructor arguments
    ---------------------
        S (spec): specification from the RWS class

    Optimal arguments
    -----------------
        none

    Attributes
    ----------
        TODO

    Methods
    -------
        TODO
    """

    def __init__(self, S, template={}):
        self.spec = S

        self.p = template.get('parameters', 0)
        self.par_len = len(self.p)
        self.par = template.get('values', np.zeros(self.par_len))
        # TODO: give a nice template (linear and quadratic),
        # if none is supplied
        self.k_sym = template.get(
            'controller', sp.zeros(len(self.spec.input), 1))
        self.V_sym = template.get('certificate', 0)

        self.dV_sym = sp.diff(self.V_sym, [self.spec.var])
        self.dtV_sym = sp.Matrix(self.dV_sym).dot(
            self.spec.system).subs(zip(self.spec.input, self.k_sym))
        self.dV_sym = None
        self.dVfun = None

        self.fitness = None

        # Create the conditions dependent on the specification
        self.conditions_sym = self.spec.create_conditions(self)
        # Create symbolic fitness function for this individual.
        self.create_fitness()

        # Create conditions for the current parameters
        self.conditions = None
        self.substitute_parameters()

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
                self.p + self.spec.extended_var, f,
                modules=[{'amax': custom_amax, 'amin': custom_amin}, "numpy"])
            result = result + (function,)
        self.fitness = result
