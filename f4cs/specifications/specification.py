# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:20:07 2021

@author: Cees F. Verdier
"""

import numpy as np
import sympy as sp
from .. import verification

# TODO: initalize with [0]*n creates a n times a copy. This can create odd
# effects


def custom_amax(x, **kwargs):
    """Alternative maximum for lambdafication."""
    return np.maximum(x[0], x[1])


def custom_amin(x, **kwargs):
    """Alternative minimum for lambdafication."""
    return np.minimum(x[0], x[1])


class Specification:
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

    def __init__(self, options, number_conditions):
        # Initialize a random generator
        self.rng = np.random.default_rng()  # Initialize random generator

        self._number_conditions = number_conditions  # number of conditions
        self.options = options
        self.var = options['variables']
        self.input = options['inputs']
        self.system = options['system']
        self.uncertainties = options.get('uncertainties', ())

        # Robustness buffer
        self.epsilon = options.get("epsilon", 0)
        # Arbitrary small constant. Default 0.01
        self.c = options.get("c", 0.01)

        # Initialize parameters with default values.
        self.number_samples = options.get("number_samples", 10000)
        self.max_samples = options.get("max_samples", 10000)

        # TODO: this list containts n copies! change
        self._data_sets = [
            np.array([[]] * self.number_samples)] * self._number_conditions
        self.condition_sets = (True,)*self._number_conditions

        # If there are uncertanties, extract the corresponding interval
        if self.uncertainties:
            self.uncertainty_list = options['uncertainty_list']
            self.extended_var = self.uncertainties + self.var
            self.add_disturbances()
        else:
            self.extended_var = self.var

        # Select the SMT solver
        self.smt_options = options.get("smt_options",
                                       {"solver": 'Z3'})
        smt_solver = self.smt_options.get("solver", 'dReal')
        smt_dictionary = {'dReal': verification.Dreal,
                          'Z3': verification.Z3}
        self.verifier = smt_dictionary[smt_solver](self.smt_options)

        # TODO: this list containts n copies! change
        self.verification_result = [
            {"sat": False, "violation": [0] * len(self.extended_var)}
        ] * self._number_conditions

    def sample_set(self, set_list, number_samples=None):
        """Sample a set defined by a set of intervals."""
        # If there is no argument of numpsamp, use atribute
        n = len(set_list)
        if number_samples is None:
            number_samples = self.number_samples
        array = np.array(set_list)
        samples = (array[:, 1] - array[:, 0]) * self.rng.random(
            (number_samples, n)
        ) + array[:, 0]
        return samples

    def sample_set_complement(self, outer_list, inner_list,
                              number_samples=None):
        """Sample a set of the form Outer not(Inner)."""
        n = len(outer_list)
        # If there is no argument of numpsamp, use atribute
        if number_samples is None:
            number_samples = self.number_samples
        outer_array = np.array(outer_list)
        inner_array = np.array(inner_list)
        samples = (outer_array[:, 1] - outer_array[:, 0]) * self.rng.random(
            (number_samples, n)
        ) + outer_array[:, 0]

        # Cut out a set
        # select a dimension
        sdims = self.rng.integers(n, size=number_samples)
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

        fit_data = [solution.fitness[i](*par, *self._data_sets[i].T)
                    for i in range(self._number_conditions)]

        norm_fit_data = np.array(
            [self.normalized_fitness(data) for data in fit_data])
        weights = self.fitness_weights(norm_fit_data)

        return np.sum(weights*norm_fit_data)/self._number_conditions

    def parameter_fitness(self, parameter, solution):
        """Given a parameter vector, compute the sample-based fitness."""
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

    # TODO:store the sample fitness of a candidate after parameter optimization
    def fitness(self, solution):
        """Compute the full fitness of a candidate solution."""
        return (self.sample_fitness(solution) + self.SMT_fitness(solution)) / 2

    def create_symbolic_interval(self, variables, set_list, open_set=False):
        """Create a symbolic interval.

        Given a set of symbolic variables and an list of interval bounds,
        create a symbolic interval. If the boolean open_set is True, it returns
        an open interval, otherwise a closed interval.

        """
        symbolic_set = sp.And()
        if not open_set:
            for i, var in enumerate(variables):
                symbolic_set = sp.And(
                    symbolic_set, sp.And(var >= set_list[i][0],
                                         var <= set_list[i][1])
                )
        else:
            for i, var in enumerate(variables):
                symbolic_set = sp.And(
                    symbolic_set, sp.And(var > set_list[i][0],
                                         var < set_list[i][1])
                )
        return symbolic_set

    def create_conditions(self, solution):
        """Create the conditions to be verified with an SMT solver."""
        return (True,)*self._number_conditions

    def add_disturbances(self):
        """Add disturbances.

        prepends specification sets with disturbances.
        """

        aux_data = self.sample_set(self.uncertainty_list)
        aux_set = self.create_symbolic_interval(self.uncertainties,
                                                self.uncertainty_list)
        # Overload system set, and sample sets
        self.add_condition_sets([aux_set]*self._number_conditions)
        self.add_data_sets([aux_data]*self._number_conditions)

    def add_condition_sets(self, condition_sets):
        """Add conditions to the specification."""
        self.condition_sets = tuple(
            [sp.And(condition, condition_sets[i])
             for i, condition in enumerate(self.condition_sets)]
        )

    def add_data_sets(self, new_data_sets):
        """Add data sets to the specification."""
        for i in range(self._number_conditions):
            self._data_sets[i] = np.append(self._data_sets[i],
                                           new_data_sets[i], axis=1)

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
                self.condition_sets[i],
                self.extended_var
            )
            # If there is a counter example, sample the data and append.
            if isinstance(self.verification_result[i]["violation"], tuple):
                self._data_sets[i] = np.append(
                    self._data_sets[i],
                    [self.verification_result[i]["violation"]],
                    axis=0,
                )
                # Saturate w.r.t. the maximum number of samples. (FIFO)
                if len(self._data_sets[i]) > self.max_samples:
                    self._data_sets[i] = self._data_sets[i][-self.max_samples:]
