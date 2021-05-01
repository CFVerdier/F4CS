# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 23:31:18 2021

@author: Cees F. Verdier
"""

from f4cs.synthesis.gggp.grammar import GrammarTree
from f4cs.synthesis.synthesis import Synthesis
from f4cs.candidates import Solution
from joblib import Parallel, delayed
import sympy as sp
import copy
import cma
import time


class GrammarGuidedGeneticProgramming(Synthesis):
    """Template-based synthesis."""

    def __init__(self, options={}):
        Synthesis.__init__(self, options)
        self.sigma0 = options.get('sigma_0', 0.5)
        self.grammar = options['grammar']
        self.number_individuals = options.get('number_individuals', 6)
        self.max_depth = options.get('max_depth', 2)

        self.number_genes = len(self.grammar.start)
        self.population = [GrammarTree(self.grammar, self.max_depth)
                           for _ in range(0, self.number_individuals)]
        self.sample_fitness = [0]*self.number_individuals
        self.candidates = [0]*self.number_individuals

    def select_operator_point(self, tree, positions):
        """Find a node to apply a genetic operator on."""
        index = tree.rng.integers(0, len(positions))
        return positions[index]

    def find_position_gene(self, tree, function, gene_int):
        """find the positions of a function on a certain gene"""
        positions = tree.get_position(function)
        return list(filter(lambda x: x[0] == gene_int, positions))

    # TODO: use a mution/crossover rate PER gene.
    def crossover(self, tree_1, tree_2):
        """Crossover operator."""
        # Find matching nonterminals
        t1 = copy.copy(tree_1)
        t2 = copy.copy(tree_2)
        nonterminals = self.grammar.nonterminals
        for gene_index in range(0, self.number_genes):
            number_nonterminals = len(nonterminals)
            positions_1 = [0]*number_nonterminals
            positions_2 = [0]*number_nonterminals
            overlap_indices = []
            for i, nonterminal in enumerate(nonterminals):
                function = sp.Function(nonterminal.capitalize())
                positions_1[i] = self.find_position_gene(t1, function,
                                                         gene_index)
                positions_2[i] = self.find_position_gene(t2, function,
                                                         gene_index)
                # Check if there is an overlap
                if all([positions_1[i], positions_2[i]]):
                    overlap_indices.append(i)
            if any(overlap_indices):
                nonterminal_choice = t1.rng.choice(overlap_indices)
                point_1 = self.select_operator_point(
                    t1, positions_1[nonterminal_choice])
                point_2 = self.select_operator_point(
                    t2, positions_2[nonterminal_choice])
                subtree = t2._positions[point_2]
                t1.replace(point_1, subtree)
        return t1

    def mutation(self, tree, max_depth):
        """Mutation operator."""
        t = copy.copy(tree)
        nonterminals = self.grammar.nonterminals
        number_nonterminals = len(nonterminals)
        for gene_index in range(0, self.number_genes):
            positions = [0]*number_nonterminals
            # Find which nonterminals are present
            present_indices = []
            for i, nonterminal in enumerate(nonterminals):
                function = sp.Function(nonterminal.capitalize())
                positions[i] = self.find_position_gene(t, function, gene_index)
                if positions[i]:
                    present_indices.append(i)
            if any(present_indices):
                nonterminal_choice = t.rng.choice(present_indices)
                point = self.select_operator_point(
                    t, positions[nonterminal_choice])
                grammar = copy.copy(self.grammar)
                grammar.start = str(nonterminals[nonterminal_choice])
                subtree = GrammarTree(grammar, max_depth)
                t.replace(point, subtree.get_tree())
        return t

    def optimize_parameters(self, template, spec):
        """Optimize parameters based on CMAES."""
        candidate = Solution(spec, template)
        # Optimize using CMAES
        res = cma.fmin(spec.parameter_fitness, candidate.par, self.sigma0,
                       args={candidate, }, options={'verbose': -9})
        return [candidate, res[0], res[1]]

    def synthesis(self, spec):
        """Synthesize controller."""
        # Create a candidate object based on the template.
        self.issolution = False

        result = Parallel(n_jobs=6)(delayed(self.optimize_parameters)
                                    (individual.create_template(), spec)
                                    for individual in self.population)
        self.candidates = [res[0] for res in result]
        # Subsitute best parameters back
        [individual.substitute_values(result[i][1]) for i, individual in
         enumerate(self.population)]

        self.sample_fitness = [res[2] for res in result]

        # for i, individual in enumerate(self.population):
        #     candidate = Solution(spec, individual.create_template())
        #     # Optimize using CMAES
        #     res = cma.fmin(spec.parameter_fitness, candidate.par,
        # self.sigma0,
        #                    args={candidate, }, options={'verbose': -1})
        #     optimal_values = res[0]
        #     self.sample_fitness[i] = res[1]
        #     self.candidates[i] = candidate
        #     # Subsitute best parameters back
        #     individual.substitute_values(optimal_values)

        # Keep iterating until a maximum of iterations is met, or a solution
        # # is found
        # self.iteration = 0
        # t_start = time.time()
        # while self.issolution is False and \
        #         (self.iteration < self.max_iterations):
        #     self.iteration += 1
        #     # Optimize parameters
        #     cma.fmin(spec.parameter_fitness, candidate.par, self.sigma0,
        #              args={candidate, }, options={'verbose': -1})
        #     # Verify candidate solution
        #     spec.verify(candidate)
        #     verification_booleans = [result['sat']
        #                              for result in spec.verification_result]
        #     self.issolution = all(verification_booleans)

        # self.synthesis_time = time.time()-t_start
        # if self.issolution is True:
        #     print('Solution found in {}'.format(self.iteration)
        #           + ' iterations')
        #     print('Synthesis time: {}'.format(self.synthesis_time)
        #           + ' seconds')
        # else:
        #     print('No solution found in {}'.format(self.iteration)
        #           + ' iterations')
        #     print('Synthesis time: {}'.format(self.synthesis_time)
        #           + ' seconds')
        # # Update log of synthesis times and iterations
        # self.iterations.append(self.iteration)
        # self.synthesis_times.append(self.synthesis_time)

        return self.population


# # Demo on grammar useage #
# # Production rules
# P = {'pol': ('pol+pol', 'const*mon'),
#      'mon': ('mon*mon', 'vars'),
#      'vars': ('x1', 'x2')
#      }
# P_terminal = {'pol': ('const*mon',),
#               'mon': ('vars',),
#               'vars': ('x1', 'x2')
#               }
# # Special production rules for constants
# P_constants = {'const': ('Real(-1,1)',)}

# start = ('pol+const', 'mon')
# grammar = Grammar(start, P, P_terminal, P_constants)
# max_depth = 3
# number_individuals = 6
# t1 = GrammarTree(grammar, max_depth)
# t2 = GrammarTree(grammar, max_depth)

# print(str(t1))
# print(str(t2))

# # Create a template that can be used in the Solution class.
# template_t1 = t1.create_template()
# # Update parameters (e.g. through optimization)
# template_t1['values'][0] = 10
# # Substitute back into tree
# t1.substitute_values(template_t1['values'])

# print(str(t1))

# # Demo for GGGP/genetic operators
# options = {'grammar': grammar,
#            'max_depth': max_depth,
#            'number_individuals': number_individuals}
# gggp = GrammarGuidedGeneticProgramming(options)

# t3 = gggp.crossover(t1, t2)
# t4 = gggp.mutation(t1, max_depth)
# print(str(t3))
# print(str(t4))
