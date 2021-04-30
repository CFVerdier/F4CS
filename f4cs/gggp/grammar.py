# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 16:13:46 2021

@author: Cees F. Verdier
"""

import sympy as sp
import numpy as np


class Grammar:
    """Grammar class."""

    def __init__(self, P, PT, start):
        self.nonterminals = tuple(P.keys())  # nonterminals
        self.P = P  # Production rules
        self.PT = PT  # terminal rules
        self.start = sp.sympify(start, evaluate=False)  # Starting tree


class GrammarTree:
    """A tree based on a grammar"""

    def __init__(self, grammar, max_depth):
        self.grammar = grammar
        self._tree = None
        self._positions = None
        self._phenotype = None
        self.rng = np.random.default_rng()  # Initialize random generator
        self.set_tree(self.grammar.start)
        self.grow_tree(max_depth)

    def __repr__(self):
        return str(self._tree)

    def __str__(self):
        return str(self.phenotype)

    def set_tree(self, new_tree):
        """Setter for the tree attribute"""
        self._tree = sp.sympify(new_tree, evaluate=False)
        self._positions = self.create_position_dictionary()
        self._phenotype = self.create_phenotype()

    def get_tree(self):
        """Getter for the tree attribute"""
        return self._tree

    def grow_tree(self, max_depth):
        """Grow a tree based on the grammar."""
        nonterminals_present = self.has_nonterminals()
        depth = 0
        # While the maximum depth is not obtained and while there are still
        # nonterminals, unfold using the default production rules
        while (depth <= max_depth) and nonterminals_present:
            self.unfold(self.grammar.P)
            nonterminals_present = self.has_nonterminals()
            depth = depth + 1
        # Maximum depth is reached: only use non-recursive rules
        while nonterminals_present:
            self.unfold(self.grammar.PT)
            nonterminals_present = self.has_nonterminals()
        return

    def unfold(self, production_rules):
        """Unfold one iteration of nonterminals"""
        nonterminals = self.grammar.nonterminals
        all_positions = [self.get_position(nonterminal)
                         for nonterminal in nonterminals]
        for i, nonterminal in enumerate(nonterminals):
            rules = production_rules[nonterminal]
            positions = all_positions[i]
            function = sp.Function('{}'.format(nonterminal.capitalize()))
            for position in positions:
                rule_index = self.rng.integers(0, high=len(rules))
                self.replace(
                    position,
                    function(sp.sympify(rules[rule_index], evaluate=False))
                )

    def create_phenotype(self):
        """Obtain the phenotype/function form."""
        phenotype = self.get_tree()
        for nonterminal in self.grammar.nonterminals:
            function = sp.Function('{}'.format(nonterminal.capitalize()))
            phenotype = phenotype.replace(function, lambda x: x)
        self.phenotype = phenotype

    def replace(self, position, new_tree):
        """
        Replace the subtree at the given position with a new subtree

        Parameters
        ----------
        position : list
            position list of the subtree to replace
        new_tree : string or sympy expression
            new tree to replace with the subtree at position
        """
        # TODO: use self.positions to create a tree? Then we can also use that
        # to easily replace subtrees
        tree = self.get_tree()
        new_tree = sp.sympify(new_tree, evaluate=False)
        depth = len(position)
        if depth > 0:
            func_list = [0]*depth
            arg_list = [0]*depth
            arg_list[0] = list(tree.args)
            func_list[0] = tree.func
            # For each layer, extract the functions and its arguments
            for i in range(1, depth):
                subtree = arg_list[i-1][position[i-1]]
                func_list[i] = subtree.func
                arg_list[i] = list(subtree.args)
            # At the desired position, replace the arguments
            arg_list[depth-1][position[depth-1]] = new_tree
            # In reverse, substitute arguments back into the functions
            for i in reversed(range(0, depth-1)):
                arg_list[i][position[i]] = func_list[i+1](*arg_list[i+1],
                                                          evaluate=False)
            # Construct the final function
            self.set_tree(func_list[0](*arg_list[0], evaluate=False))

        else:
            self.set_tree(new_tree)

    def create_subtree_dictionary(self, tree, position_dict={}, position=()):
        """ Get all positions of the elements in the set elements"""
        position_dict[position] = tree
        if tree.args != []:
            for pos, arg in enumerate(tree.args):
                subposition = position+(pos,)
                self.create_subtree_dictionary(arg, position_dict, subposition)
        return position_dict

    def create_position_dictionary(self):
        """Wrapper function for get_subtree_dictionary for the tree"""
        tree = self.get_tree()
        self.positions = self.create_subtree_dictionary(tree)

    def get_position(self, root):
        """Get the positions in which the root occur.

        The function distinguishes symbols and functions. Examples:
        tree.get_position('pol')
        tree.get_position(sp.Add)
        """
        root = sp.sympify(root, evaluate=False)
        if root.is_Symbol:
            positions = [position for position, value in self.positions.items()
                         if str(value) == str(root)]
        else:
            positions = [position for position, value in self.positions.items()
                         if value.func == root]
        return positions

    def has_nonterminals(self):
        """check whether the tree still has nonterminals"""
        return any([self.get_position(nonterminal) for nonterminal in
                    self.grammar.nonterminals])


# Production rules
P = {'pol': ('pol+pol', 'mon'),
     'mon': ('mon*mon', 'vars'),
     'vars': ('x1', 'x2')
     }
P_terminal = {'pol': ('mon',),
              'mon': ('vars',),
              'vars': ('x1', 'x2')
              }

start = 'pol'
grammar = Grammar(P, P_terminal, start)
max_depth = 3
t1 = GrammarTree(grammar, 3)

print(str(t1))
