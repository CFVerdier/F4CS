# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 17:24:46 2018

@author: ceesv
"""

import sympy as sp
import numpy as np
from sympy.parsing.sympy_parser import parse_expr


class Tree(object):
    """Represent a Tree structure."""

    def __init__(self, tp=0, chlds=None):
        self.type = tp
        self.children = chlds
#        if chlds is not None:
#            self.children = {x: chlds[x] for x in list(chlds.keys()) }

    def __str__(self):
        """Convert object to string."""
        childstring = ''
        if self.children is not None:
            for child in self.children:
                childstring = childstring + str(child) + ','
            childstring = childstring.rstrip(',')
        return str(self.type) + "(" + childstring + ")"

    def depth(self):
        """Return tree depth."""
        childdepth = 0
        if self.children is not None:
            for child in self.children:
                childdepth = np.max([child.depth(), childdepth])
        return 1 + childdepth

    def evaluate(self):
        """Evaluate tree."""
        string = str(self.type) + '('
        if self.children is not None:
            for child in self.children:
                string = string + child.evaluate() + ','
            string = string.rstrip(',')
        return string + ')'

    def functionForm(self):
        """Give the string form."""
        string = ''
        if self.children is None:
            string = string + str(self.type)
            string = string.lstrip('(')
            string = string.rstrip(')')

        else:
            for child in self.children:
                string = string + str(child.functionForm()) + self.type
            string = string.rstrip(self.type)
        return parse_expr(string)


class GTree(object):
    """Represent a Grammar Tree structure."""

    def __init__(self, nonterminal, tree):
        self.nonterminal = nonterminal
        self.tree = tree


class GPfunctionals(object):
    def __init__(self):
        return

    def plus(self, args):
        children = {}
        for arg in range(1, len(args)+1):
            children[arg] = Tree(args[arg-1])
        return Tree('+', children)

    def times(self, args):
        children = {}
        for arg in range(1, len(args)+1):
            children[arg] = Tree(args[arg-1])
        return Tree('*', children)

    def f(self, name, args=None):
        """ Generic gp function"""
        children = {}
        for arg in range(1, len(args)+1):
            children[arg] = Tree(args[arg-1])
        return Tree(name, children)


class Grammar(object):
    """Grammar class."""

    def __init__(self, P, PT, S):
        self.n = P.keys
        self.P = {}
        self.PT = {}
        self.create_internal_rules(P, PT)
        self.s = GTree(S, S)

    def create_internal_rules(self, P_user, PT_user):
        for n in P_user:
            self.P[n] = [GTree(n, entry) for entry in P_user[n]]
        for n in PT_user:
            self.PT[n] = [GTree(n, entry) for entry in PT_user[n]]


# Demo tree
x, y, z = sp.symbols("x y z")

t1 = Tree('+', [Tree(x), Tree(y)])
t2 = Tree('*', [Tree('1'), Tree('2'), Tree('+', [Tree('3'), Tree('4')])])


print(t1.evaluate())
print(t1.functionForm())
#    print(t2)
print(t2.evaluate())

# use function
print(t1.functionForm().subs([(x, 2), (y, 1)]))

x1, x2 = sp.symbols("x1 x2")

# P = {'expr':['pol'],
#     'pol':[gp.plus(['pol','pol']),'pol'],
#     'mon':[gp.times(['var','mon']),gp.times(['mon','mon'])],
#     'var':[x1,x2]}
P = {'pol': [Tree('+', ['pol', 'pol']), 'mon'],
     'mon': [Tree('*', ['mon', 'mon']), 'var'],
     'var': [x1, x2]}
PT = {'pol': ['pol'],
      'mon': ['var'],
      'var': [x1, x2]}

start = Tree('pol')

# create grammar
g1 = Grammar(P, PT, start)

# start tree
t1 = g1.s


# Note that we need to copy this object (works as a pointer)
# t1 = copy.deepcopy(start)
#
# start = Tree('pol')

# Give a list of nonterminals, and objects are created automatically.
