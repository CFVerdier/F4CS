# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:18:27 2021

@author: Cees F. Verdier
"""

from .specification import Specification
import sympy as sp

class RWS(Specification):
    """Reach-while-stay specification.

    Represents the RWS specification and implements the corresponding
    fitness function and verification. Subclass of Specification.

    Arguments
    ---------
        variables: tuple of symbolic state variables.
            Example: sympy.symbols('x1, x2')
        inputs: tuple of symbolic input variables.
            Example: sympy.symbols('u1,')
        f_sym (sympy expression): Symbolic expression of the system dynamics
        options (dictionary): dictionary with all the settings relating to the
          specification.

    Required options
    ----------------
        Slist, Ilist, Olist: a list of the lower and upperbounds of the
          safe set, initial set and goal set
        path: path where the SMT files are stored.

    Optional
    --------
        number_samples:  Number of samples. Default 100
        rdelta:         Inflation of the boundary. Default: 0.01
        gamma:          (Arbitrary) decrease of the LF. Default :0.01,
        c:              (Arbitrary) nonnegative parameter (see manual).
                          Default: 0.01
        dprecision:     Precision of dReal. Default: 0.01


    """

    def __init__(self, variables, inputs, f_sym, options):
        # Call the __init__ function of the Spec parent class first.
        Specification.__init__(self, variables, inputs, f_sym, options)

        self._number_conditions = 3  # number of RWS conditions

        S_list = self.options["Slist"]
        I_list = self.options["Ilist"]
        O_list = self.options["Olist"]
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
        S_not_O_data = self.sample_set_complement(S_list, O_list)
        self.data_sets = [I_data, dS_data, S_not_O_data]

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

        self.condition_set = (I_set, closed_R_not_S_set, S_not_O_set)
        # TODO: this list containts n copies! change
        self.verification_result = [None] * self._number_conditions

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        # create the conditions to verify
        con1 = solution.V_sym <= 0
        con2 = solution.V_sym > 0
        con3 = sp.Or(solution.V_sym > 0, solution.dtV_sym <= -self.gamma)

        return (con1, con2, con3)