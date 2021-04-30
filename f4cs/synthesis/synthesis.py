# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:54:28 2021

@author: Cees F. Verdier
"""


class Synthesis:
    """Main synthesis class."""

    def __init__(self, options={}):
        self.issolution = False
        self.mode = options.get('mode', 'template-based')
        self.max_iterations = options.get('max_iterations', 10)
        self.iteration = None
        self.synthesis_time = None
        # Store past number of iterations and synthesis time. Useful for
        # logging past experiments
        self.iterations = []
        self.synthesis_times = []
