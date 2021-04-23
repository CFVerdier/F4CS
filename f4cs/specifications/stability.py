# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 21:35:41 2021

@author: ceesv
"""

import numpy as np
import sympy as sp
from .specifications import Specification

print("Warning: this is an experimental module")

class Stability(Specification):
    """Stability specification.

    For polynomial systems, this specification checks whether for a domain D,
    the standard LF conditions hold.

    For nonpolynomial systems, see local_stability
    """

    def __init__(self, options):
        # Call the __init__ function of the Spec parent class first.
        number_conditions = 2  # number of local stability conditions
        Specification.__init__(self, options, number_conditions)

        D_list = self.options["Dlist"]
        self.x0 = self.options.get('x0', sp.zeros(self.n, 1))
        # decrease of the LBF. Default 0.01
        self.gamma = self.options.get("gamma", 0.01)
        self.positive_definite = self.options.get(
            "PD", self.gamma*np.sum([xi**2 for xi in self.var]))

        # Create sample sets
        D_data = self.sample_set(D_list)
        self.add_data_sets([D_data, D_data])

        # Create symbolic domains for SMT solver
        D_set = self.create_symbolic_interval(self.var, D_list)

        self.add_condition_sets((D_set, D_set))
        self.conditions = None
        self.verification_result = [None] * self._number_conditions

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        con1 = solution.V_sym >= self.positive_definite
        con2 = solution.dtV_sym <= -self.positive_definite

        return (con1, con2,)
