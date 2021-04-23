# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:54:28 2021

@author: Cees F. Verdier
"""

import cma
from .candidates import Solution


class Synthesis:
    """Main synthesis class."""

    def __init__(self, options={}):
        self.issolution = False
        self.mode = options.get('mode', 'template-based')
        self.max_iterations = options.get('max_iterations', 10)
        self.iteration = 0


class TBS(Synthesis):
    """Template-based synthesis."""

    def __init__(self, options={}):
        Synthesis.__init__(self, options)
        self.sigma0 = options.get('sigma_0', 0.5)

    def synthesis(self, spec, template={}):
        """Synthesize controller."""
        # Create a candidate object based on the template.
        candidate = Solution(spec, template)

        # Keep iterating until a maximum of iterations is met, or a solution
        # is found
        while self.issolution is False and \
                (self.iteration < self.max_iterations):
            self.iteration += 1
            # Optimize parameters
            cma.fmin(spec.parameter_fitness, candidate.par, self.sigma0,
                     args={candidate, }, options={'verbose': -1})
            # Verify candidate solution
            spec.verify(candidate)
            verification_booleans = [result['sat']
                                     for result in spec.verification_result]
            self.issolution = all(verification_booleans)

        if self.issolution is True:
            print('Solution found in {}'.format(self.iteration)
                  + ' iterations')
        else:
            print('No solution found in {}'.format(self.iteration)
                  + ' iterations')

        return candidate
