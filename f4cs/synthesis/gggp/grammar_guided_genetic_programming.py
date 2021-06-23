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
import numpy as np
import copy
import cma
import time


class GrammarGuidedGeneticProgramming(Synthesis):
    """Template-based synthesis."""

    def __init__(self, options={}):
        Synthesis.__init__(self, options)
        self.sigma0 = options.get('sigma_0', 0.5)
        self.grammar = options['grammar']
        self.number_individuals = options.get('number_individuals', 12)
        self.max_depth = options.get('max_depth', 2)

        self.tournament_size = options.get('tournament_size',
                                           int(np.min(
                                               [2, np.floor(
                                                   self.number_individuals*0.2
                                               )])))
        crossover_rate = options.get('crossover_rate', 0.8)
        mutation_rate = options.get('mutation_rate', 0.3)
        # Genetic operator rates
        self.operator_rates = {'crossover': crossover_rate,
                               'mutation': mutation_rate}
        self.cores = options.get('cpu_cores', -2)  # default is all but 1

        self.number_genes = len(self.grammar.start)
        self.population = [GrammarTree(self.grammar, self.max_depth)
                           for _ in range(0, self.number_individuals)]
        self.sample_fitness = [0]*self.number_individuals
        self.SMT_fitness = [0]*self.number_individuals
        self.candidates = [0]*self.number_individuals

        self.best = None
        self.best_fitness = None

        self.rng = np.random.default_rng()  # Initialize random generator
        self.log = {'best_fitness': [],
                    'average_fitness': []}

    def select_operator_point(self, tree, positions):
        """Find a node to apply a genetic operator on."""
        index = self.rng.integers(0, len(positions))
        return positions[index]

    def find_indices(self, a_list, element):
        """Given a list and an element, find the positions."""
        return [i for i, e in enumerate(a_list) if e == element]

    def find_position_gene(self, tree, function, gene_int):
        """find the positions of a function on a certain gene"""
        positions = tree.get_position(function)
        return list(filter(lambda x: x[0] == gene_int, positions))

    def tournament_selection(self):
        """Tournament selection"""
        selected_indices = self.rng.integers(0, high=self.number_individuals,
                                             size=self.tournament_size)
        selected = [self.population[i] for i in selected_indices]
        fitness_selected = [self.full_fitness[i] for i in selected_indices]
        best = selected[fitness_selected.index(max(fitness_selected))]
        return best

    def create_new_individual(self):
        """Create a new individual, based on selection and genetic operators.

        """
        operator_probabilities = self.rng.random(2)
        tree_1 = self.tournament_selection()
        tree_2 = self.tournament_selection()
        new_tree = tree_1
        # Apply crossover
        if operator_probabilities[0] <= self.operator_rates['crossover']:
            new_tree = self.crossover(tree_1, tree_2)
        # Apply mutation
        if operator_probabilities[1] <= self.operator_rates['mutation']:
            new_tree = self.mutation(new_tree, self.max_depth)
        return new_tree

    def create_new_population(self):
        """Create new population, based on elitism, selection, and genetic
        operators.

        """
        # Eliteism
        new_population = [0]*self.number_individuals
        new_population[0] = self.best
        # Creating a new population
        for i in range(1, self.number_individuals):
            new_population[i] = self.create_new_individual()
        self.population = new_population

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
                nonterminal_choice = self.rng.choice(overlap_indices)
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
                nonterminal_choice = self.rng.choice(present_indices)
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
        res = cma.fmin(spec.parameter_fitness,
                       candidate.par,
                       self.sigma0,
                       args={candidate, },
                       options={'verbose': -9,
                                'ftarget': -1,
                                # 'maxiter': 30,
                                'CMA_diagonal': True,
                                'timeout': 1})
        return [candidate, res[0], res[1]]

    def synthesis(self, spec):
        """Synthesize controller."""
        # Create a candidate object based on the template.
        self.issolution = False
        self.iteration = 0
        t_start = time.perf_counter()

        while ((self.issolution is False) and
                (self.iteration < self.max_iterations)):

            self.iteration += 1
            # optimize parameters
            t_start_optimization = time.perf_counter()
            result = Parallel(n_jobs=self.cores, prefer="threads")(
                delayed(self.optimize_parameters)
                (individual.create_template(), spec)
                for individual in self.population)
            self.sample_fitness = [-res[2] for res in result]
            # Subsitute best parameters back
            [individual.substitute_values(result[i][1]) for i, individual in
             enumerate(self.population)]
            t_end_optimization = time.perf_counter()-t_start_optimization
            print('Parameters optimized in: ' + str(t_end_optimization) + ' s')
            
            t_start_verify = time.perf_counter()
            # Verify candidate solution (only those with a maximum samp. fit)
            to_verify = self.find_indices(self.sample_fitness, 1.0)
            self.smt_fitness = [0]*self.number_individuals
            for i in to_verify:
                self.smt_fitness[i] = spec.SMT_fitness(self.candidates[i])

            self.full_fitness = list(
                np.add(self.sample_fitness, self.smt_fitness))
            t_end_verify = time.perf_counter()-t_start_verify
            print('Verified in: ' + str(t_end_verify) + ' s')

            best_index = self.full_fitness.index(max(self.full_fitness))
            self.best_fitness = self.full_fitness[best_index]
            self.best = self.population[best_index]

            # Log results
            self.log['best_fitness'].append(self.best_fitness)
            self.log['average_fitness'].append(np.mean(self.full_fitness))

            print('Best fitness: ' + str(self.best_fitness))

            # Check if there is a solution
            if (self.best_fitness == 2):
                self.issolution = True
                break

            # Create new population based on genetic operators.
            t_start_creation = time.perf_counter()
            self.create_new_population()
            t_end_creation = time.perf_counter()-t_start_creation
            print('New population created in: ' +
                  str(t_end_creation) + ' s')

        # After loop, check if the synthesis was succesfull
        self.synthesis_time = time.perf_counter()-t_start
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

        return self.best
