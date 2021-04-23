# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 21:35:41 2021

@author: ceesv
"""

import numpy as np
import sympy as sp
from .specifications import Specification

print("Warning: this is an experimental module")


class PracticalStability(Specification):
    """Practical Stability specifcation.

    This specification checks whether for a domain D\0,
    the standard LF conditions hold, where 0 is a set around the equilibrium.
    """

    def __init__(self, options):
        # Call the __init__ function of the Spec parent class first.
        number_conditions = 2  # number of local stability conditions
        Specification.__init__(self, options, number_conditions)

        D_list = self.options["Dlist"]
        O_list = self.options["Olist"]
        # decrease of the LBF. Default 0.01
        self.gamma = self.options.get("gamma", 0.01)

        # Create sample sets
        D_not_O_data = self.sample_set_complement(D_list, O_list)
        self.add_data_sets([D_not_O_data, D_not_O_data])

        # Create symbolic domains for SMT solver
        D_set = self.create_symbolic_interval(self.var, D_list)
        O_set = self.create_symbolic_interval(self.var, O_list)

        D_not_O_set = sp.And(D_set, sp.Not(O_set))

        self.add_condition_sets((D_not_O_set, D_not_O_set))
        self.conditions = None
        self.verification_result = [None] * self._number_conditions

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        con1 = solution.V_sym > 0
        con2 = solution.dtV_sym <= -self.gamma

        return (con1, con2,)
