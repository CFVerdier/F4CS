# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:20:07 2021

@author: ceesv
"""

import numpy as np
import sympy as sp
import verification
# import verification.dReal_communication as dReal

# TODO: initalize with [0]*n creates a n times a copy. This can create odd
# effects
# TODO: create symbolic fitness with parameters not substituted in.


def custom_amax(x, **kwargs):
    """Alternative maximum for lambdafication."""
    return np.maximum(x[0], x[1])


def custom_amin(x, **kwargs):
    """Alternative minimum for lambdafication."""
    return np.minimum(x[0], x[1])


class Spec:
    """Base class of a specification.

    Constructor arguments
    ---------------------
        variables: tuple of symbolic state variables.
            Example: sympy.symbols('x1, x2')
        inputs: tuple of symbolic input variables.
            Example: sympy.symbols('u1')
        options (dictionary): dictionary with all the specification settings.

    Attributes
    ----------
        var (list): User defined state variables
        input (list): User defined system inputs
        n (int): number system variables
        m (int): number of input variables


    Methods
    -------
        sample_set(interval, numper of samples): sample a number of random
          samples from a hyperrectangle, spanned by the interval. Returns an
          array of float64.
        sample_set_complement(OuterList,InnerList,number_samples): Samples a
          set Outerset not(Innerset), where number_samples denotes the number
          of samples, Outerlist and Innerlist the interval of the Outerset and
          Innerset respectively.

    Child classes:
        RWS: Reach while stay specification for continuous time systems

    """

    def __init__(self, variables, inputs, f_sym, options):

        self._number_conditions = 1  # number of conditions
        self.var = variables
        self.input = inputs
        self.options = options
        # Initialize parameters with default values.
        self.dprecision = self.options.get("dprecision", 0.01)
        self.number_samples = self.options.get("number_samples", 100)
        self.max_samples = self.options.get("max_samp", 1000)
        # Robustness buffer
        self.epsilon = self.options.get("epsilon", 0.1)
        # Arbitrary small constant. Default 0.01
        self.c = self.options.get("c", 0.01)
        self.smt_options = self.options.get("smt_options",
                                            {"solver": 'dReal'})

        # Select the SMT solver
        smt_solver = self.smt_options.get("solver", 'dReal')
        smt_dictionary = {'dReal': verification.Dreal,
                          'Z3': verification.Z3}
        self.verifier = smt_dictionary[smt_solver](self.smt_options)

        # Make functions of the dynamics
        self.f_sym = f_sym

        self.n = len(self.var)
        self.m = len(self.input)

        self.condition_set = (True,)

        # TODO: this list containts n copies! change
        self.verification_result = [
            {"sat": False, "violation": [0] * self.n}
        ] * self._number_conditions

        self.rng = np.random.default_rng()  # Initialize random generator

        # TODO: this list containts n copies! change
        self.data_sets = [
            np.array([[0] * self.n] * self.number_samples)
        ] * self._number_conditions

    def sample_set(self, set_list, number_samples=None):
        """Sample a set defined by a set of interfals."""
        # If there is no argument of numpsamp, use atribute
        if number_samples is None:
            number_samples = self.number_samples
        array = np.array(set_list)
        samples = (array[:, 1] - array[:, 0]) * self.rng.random(
            (number_samples, self.n)
        ) + array[:, 0]
        return samples

    def sample_set_complement(self, outer_list, inner_list,
                              number_samples=None):
        """Sample a set of the form Outer not(Inner)."""
        # If there is no argument of numpsamp, use atribute
        if number_samples is None:
            number_samples = self.number_samples
        outer_array = np.array(outer_list)
        inner_array = np.array(inner_list)
        samples = (outer_array[:, 1] - outer_array[:, 0]) * self.rng.random(
            (number_samples, self.n)
        ) + outer_array[:, 0]

        # Cut out a set
        # select a dimension
        sdims = self.rng.integers(self.n, size=number_samples)
        # percentage
        perc = 2 * self.rng.random(number_samples) - 1
        for i in range(number_samples):
            sdim = sdims[i]
            if perc[i] >= 0:
                samples[i, sdim] = inner_array[sdim, 1] + perc[i] * (
                    outer_array[sdim, 1] - inner_array[sdim, 1]
                )
            else:
                samples[i, sdim] = inner_array[sdim, 0] - perc[i] * (
                    outer_array[sdim, 0] - inner_array[sdim, 0]
                )
        return samples

    def fitness_weights(self, data):
        """Create the fitness weights."""
        weights = np.array([1])
        for i in range(1, self._number_conditions):
            weights = np.append(weights, np.floor(weights[i-1] * data[i-1]))
        return weights

    def normalized_fitness(self, data):
        """Create normalized fitness."""
        return 1 / (1 + np.linalg.norm(data))

    def sample_fitness(self, solution):
        """Compute the sample-based fitness."""
        par = solution.par

        fit_data = [solution.fitness[i](*par, *self.data_sets[i].T)
                    for i in range(self._number_conditions)]

        norm_fit_data = np.array(
            [self.normalized_fitness(data) for data in fit_data])
        weights = self.fitness_weights(norm_fit_data)

        return np.sum(weights*norm_fit_data)/self._number_conditions

    def parameter_fitness(self, parameter, solution):
        """Given a parameter vector, compute the sample-based fitness."""
        # solution.substitute_parameters(parameter)
        # Substitute parameters
        solution.par = parameter
        return -1 * self.sample_fitness(solution)

    def SMT_fitness(self, solution):
        """Compute the SMT-based fitness."""
        # Call SMT solver to verify the conditions.
        self.verify(solution)
        # Total SMT fitness is equal to the number of true statements,
        # normalized to number of conditions
        return (
            np.sum([int(n["sat"]) for n in self.verification_result])
        ) / self._number_conditions

    # TODO: if we store the sample fitness of a candidate after parameter
    # optimization
    def fitness(self, solution):
        """Compute the full fitness of a candidate solution."""
        return (self.sample_fitness(solution) + self.SMT_fitness(solution)) / 2

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        return (True,)

    def verify(self, solution):
        """Verify the specification using an SMT solver."""
        # Create the conditions to verify
        solution.substitute_parameters()
        # For each condition, verify the condition with dReal
        for i in range(0, self._number_conditions):
            # TODO: clean up passing a file name. Now this is only relevant for
            # dReal, but is passed for any SMT solver.
            if hasattr(self.verifier, 'file_name'):
                self.verifier.file_name = "con{}".format(i + 1)
            self.verification_result[i] = self.verifier.verify(
                solution.conditions[i],
                self.condition_set[i],
                self.var
            )
            # If there is a counter example, sample the data and append.
            if isinstance(self.verification_result[i]["violation"], tuple):
                self.data_sets[i] = np.append(
                    self.data_sets[i],
                    [self.verification_result[i]["violation"]],
                    axis=0,
                )
                # Saturate w.r.t. the maximum number of samples. (FIFO)
                if len(self.data_sets[i]) > self.max_samples:
                    self.data_sets[i] = self.data_sets[i][-self.max_samples:]


class RWS(Spec):
    """Reach-while-stay specification.

    Represents the RWS specification and implements the corresponding
    fitness function and verification. Subclass of Spec.

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
        Spec.__init__(self, variables, inputs, f_sym, options)

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


# var_list = x1, x2 = sp.symbols('x1,x2')
# input_list = u1, = sp.symbols('u1,')
# f = sp.Or(x1+x2*u1 <= 2, x1 > 0)
# c = 0.001
