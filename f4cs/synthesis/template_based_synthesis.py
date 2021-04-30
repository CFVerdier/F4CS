# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 23:07:46 2021

@author: Cees F. Verdier
"""

from .synthesis import Synthesis
from ..candidates import Solution
import cma
import time


class TemplateBasedSynthesis(Synthesis):
    """Template-based synthesis."""

    def __init__(self, options={}):
        Synthesis.__init__(self, options)
        self.sigma0 = options.get('sigma_0', 0.5)

    def synthesis(self, spec, template={}):
        """Synthesize controller."""
        # Create a candidate object based on the template.
        self.issolution = False
        candidate = Solution(spec, template)

        # Keep iterating until a maximum of iterations is met, or a solution
        # is found
        self.iteration = 0
        t_start = time.time()
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

        self.synthesis_time = time.time()-t_start
        if self.issolution is True:
            print('Solution found in {}'.format(self.iteration)
                  + ' iterations')
            print('Synthesis time: {}'.format(self.synthesis_time)
                  + ' seconds')
        else:
            print('No solution found in {}'.format(self.iteration)
                  + ' iterations')
            print('Synthesis time: {}'.format(self.synthesis_time)
                  + ' seconds')
        # Update log of synthesis times and iterations
        self.iterations.append(self.iteration)
        self.synthesis_times.append(self.synthesis_time)

        return candidate
