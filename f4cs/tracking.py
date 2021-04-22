# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 21:35:41 2021

@author: ceesv
"""

import numpy as np
import sympy as sp
from .specifications import Spec

# TODO: update to new specification format


class Tracking(Spec):
    """Reach-while-stay specification of a tracking controller.

    Implements the RWS variant of the tracking controller.
    """

    def __init__(self, variables, inputs, f_sym, options, auxillary_variables):
        # Call the __init__ function of the Spec parent class first.
        Spec.__init__(self, variables, inputs, f_sym, options)

        self.aux_var = auxillary_variables
        self._number_conditions = 4  # number of RWS conditions

        S_list = self.options["Slist"]
        I_list = self.options["Ilist"]
        O_list = self.options["Olist"]
        T_list = self.options["Tlist"]

        self.t_max = T_list[0][1]
        self.alpha = self.options.get("alpha", 1)
        self.beta = self.options.get("beta", 0.05)
        self.gamma = - np.log(self.beta/self.alpha)/self.t_max

        # decrease of the LBF. Default 0.01
        self.gamma = self.options.get("gamma", 0.01)

        # Create an inflated safe set to create a conservative boundary set
        r_delta = self.options.get("rdelta", 0.01)  # Default =0.01
        R_list = [
            [S_list[i][0] - r_delta, S_list[i][1] + r_delta]
            for i in range(0, self.n)
        ]

        # Create sample sets
        I_data = self.sample_set(I_list)
        dS_data = self.sample_set_complement(R_list, S_list)
        O_data = self.sample_set(O_list)
        S_not_O_data = self.sample_set_complement(S_list, O_list)
        self.data_sets = [I_data, dS_data, O_data, S_not_O_data]

        # Create symbolic domains for SMT solver
        S_set = sp.And()
        R_set = sp.And()
        I_set = sp.And()
        O_set = sp.And()
        for i in range(0, self.n):
            S_set = sp.And(
                S_set, sp.And(self.var[i] >= S_list[i][0],
                              self.var[i] <= S_list[i][1])
            )
            I_set = sp.And(
                I_set, sp.And(self.var[i] >= I_list[i][0],
                              self.var[i] <= I_list[i][1])
            )
            O_set = sp.And(
                O_set, sp.And(self.var[i] >= O_list[i][0],
                              self.var[i] <= O_list[i][1])
            )
            S_not_O_set = sp.And(S_set, sp.Not(O_set))
            # Create the closure of R\S to create a conservative boundary of S.
            R_set = sp.And(
                R_set, sp.And(self.var[i] >= R_list[i][0],
                              self.var[i] <= R_list[i][1])
            )
            S_open_set = sp.And(
                S_set, sp.And(self.var[i] > S_list[i][0],
                              self.var[i] < S_list[i][1])
            )
            closed_R_not_S_set = sp.And(R_set, sp.Not(S_open_set))

        self.condition_set = (I_set, closed_R_not_S_set, O_set, S_not_O_set)
        # TODO: this list containts n copies! change
        self.verification_result = [None] * self._number_conditions

        # TODO: hackish: Overload the system variables and condition sets to
        # include the auxillary variables
        T_data = self.sample_set(T_list)
        T_set = sp.And()
        for i in range(0, len(self.aux_var)):
            T_set = sp.And(
                T_set, sp.And(self.aux_var[i] >= T_list[i][0],
                              self.aux_var[i] <= T_list[i][1])
            )

        new_conditions = [sp.And(condition, T_set)
                          for condition in self.condition_set]
        # Overload system set, variables, and sample sets
        self.condition_set = tuple(new_conditions)
        self.var = self.var + self.aux_var
        for i in range(self._number_conditions):
            self.data_sets[i] = np.append(self.data_sets[i],
                                          T_data, axis=1)

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        con1 = solution.V_sym <= self.alpha
        con2 = solution.V_sym > self.alpha
        con3 = solution.V_sym <= self.beta
        con4 = sp.Or(solution.V_sym > self.alpha,
                     solution.dtV_sym <= -0.001)

        return (con1, con2, con3, con4)
