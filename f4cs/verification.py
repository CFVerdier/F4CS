# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 13:46:24 2021

@author: ceesv
"""

import numpy as np
import collections
import os.path
import platform
import z3

if platform.system() == 'Windows':
    from os_support._windows_support import call_dReal
elif (platform.system() == 'Linux') or (platform.system() == 'Darwin'):
    from ._unix_support import call_dReal
else:
    raise ImportError("The verification method does not support this OS")


class Verification:
    """SMT-based verification of certificate functions.

    Solvers: dReal, Z3.
    """

    def __init__(self, options={}):
        self.options = options

        self.sym_dict = collections.OrderedDict({"Add": "+",
                                                 "Mul": "*",
                                                 "Pow": "^",
                                                 "StrictGreaterThan": ">",
                                                 "GreaterThan": ">=",
                                                 "StrictLessThan": "<",
                                                 "LessThan": "<=",
                                                 "Implies": "=>",
                                                 "And": "and",
                                                 "Or": "or",
                                                 "Not": "not",
                                                 "Max": "max",
                                                 "Abs": "abs",
                                                 "Unequality": "!",
                                                 "Equality": "="})

        self.sym_dict_exceptions = collections.OrderedDict(
            {"Max": "max", "Min": "min"})
        self.num_dict = ["Float", "Integer", "Zero", "One",
                         "Symbol", "NegativeOne", "Rational", "Half"]

    def symbolic_name_to_lisp(self, fun):
        """Translate function between symbolic python and SMT2 function names.

        Part of the symbolic_to_lisp translator.
        """
        expr = fun
        for word in self.sym_dict:
            expr = expr.replace(word, self.sym_dict[word])
        return expr

    def symbolic_to_lisp(self, expr):
        """Translate function from symbolic python to SMT2 expressions."""
        name = expr.func.__name__
        if name in self.num_dict:
            sform = str(expr)
        elif name in self.sym_dict_exceptions:
            sform = "(" + self.symbolic_name_to_lisp(name)
            sform = sform + " " + self.ymbolic_to_lisp(expr.args[0])
            for arg in expr.args[1:-1]:
                sform = sform + " (" + self.symbolic_name_to_lisp(name) + \
                    " " + self.symbolic_to_lisp(arg)
            sform = sform + " " + self.symbolic_to_lisp(expr.args[-1])
            for arg in expr.args[1:]:
                sform = sform + ")"
        else:
            sform = "(" + self.symbolic_name_to_lisp(name)
            for arg in expr.args:
                sform = sform + " " + self.symbolic_to_lisp(arg)
            sform = sform + ")"
        return sform

    def make_SMT2_string(self, symbolic_expr, domain, symbol_list):
        """Write an SMT file in string format.

        The SMT file is used to check whether inequality symbolic_expr is
        satisfied using dReal.
        """
        # Write settings in a string
        string = ""
        if self.solver == 'dReal':
            string = "(set-logic QF_NRA)\n"
            string = string + "(set-info :precision " + \
                str(self.dreal_precision)+" )\n"
        # add declaration of constants
        for var in symbol_list:
            string = string + "(declare-fun " + str(var) + " () Real)\n"
        # define domain
        string = string + "(assert " + self.symbolic_to_lisp(domain) + ")\n"
        # Add inequality
        string = string + \
            "(assert (not  " + self.symbolic_to_lisp(symbolic_expr) + "))\n"
        # check satisfiability
        string = string + "(check-sat)\n"
        string = string + "(exit)"
        return string

class OptionsError(Exception):
    pass

class Dreal(Verification):
    """dReal based verification.

    Class the SMT solver dReal

    Options:
        dreal_precision: precision used in dReal. Default:"0.001"
        t_max: maximum time for dReal to run. Terminated if surpassed. Default:
            None.
    """

    def __init__(self, options={}):
        Verification.__init__(self, options)
        self.solver = "dReal"
        self.path = self.options['path']
        self.dreal_precision = self.options.get("dreal_precision", 0.001)
        self.t_max = self.options.get("dreal_precision", None)
        self.file_name = self.options.get("file_name", "file.smt2")

        # Check for unix-based OSes whether the dReal path is supplied
        if (platform.system() == 'Linux') or (platform.system() == 'Darwin'):
            try:
                self.dReal_path = self.options["dReal_path"]
            except KeyError:
                raise OptionsError('Please supply the dReal path to ' +
                                      'options under the key "dReal_path"')

    def call(self):
        """Call dReal.

        Returns a dictionary containing whether the inequality is satisfied,
        if not, the violations stored as a string, and a time-out flag.
        """
        # Check file suffix
        if not self.file_name.endswith('.smt2'):
            self.file_name = self.file_name + '.smt2'

        # Initialize results
        result = {}
        print("Calling dReal")
        try:
            # dReal call, OS specific
            outputdReal = call_dReal(self)
        except Exception:
            outputdReal = 'time-out'
            result['time-out'] = True
            print("dReal time-out ({} seconds"
                  " or other unexpected result.)".format(self.t_max))
        # Process data: unsat = 1, delta-sat = 0 (unsat proofs the inequality)
        if outputdReal == 'time-out':
            result['sat'] = False
            result['violation'] = None
        elif outputdReal == 'unsat\n':
            result['sat'] = True
            result['time-out'] = False
            result['violation'] = None
            print('Inequality Satisfied')
        else:
            result['sat'] = False
            result['time-out'] = False
            result['violation'] = outputdReal
            print('Inequality not verified')
        return result

    def read_violations(self, dreal_result, symbol_list):
        """Read violations for the output.

        Note that accuracy might be lost due to the use of eval()
        """
        string = dreal_result
        for var in symbol_list:
            string = string.replace(str(var), '"'+str(var)+'"')
        array = string.splitlines()
        array = array[1:]
        # Store point as a dictionary
        viol_dict = {}
        for entry in array:
            d = eval('{'+entry + '}')
            viol_dict.update(d)
        # Take the first vertex as counter example
        new_sample = ()
        for var in symbol_list:
            #        new_sample = new_sample +(viol_dict[str(var)][0],)
            new_sample = new_sample + (np.mean(viol_dict[str(var)]),)
        return new_sample

    def verify(self, symbolic_expr, domain, symbol_list):
        """Verify using dReal.

        Main wrapper function of the verification class.
        """
        # Check file suffix
        if not self.file_name.endswith('.smt2'):
            self.file_name = self.file_name + '.smt2'
        # Write SMT file
        string = self.make_SMT2_string(symbolic_expr,
                                       domain, symbol_list)
        with open(os.path.join(self.path, self.file_name), 'w+') as f:
            f.write(string)

        print("SMT2 File exported at "
              + os.path.join(self.path, self.file_name))

        # Call dReal
        result = self.call()
        # Read violations
        if result['violation'] is not None and result['time-out'] is False:
            result['violation'] = self.read_violations(
                result['violation'], symbol_list)
            print("Violation at" + str(result['violation']))
        return result


class Z3(Verification):
    """Z3 based verification.

    Class the SMT solver Z3
    """

    def __init__(self, options={}):
        Verification.__init__(self, options)
        self.solver = "Z3"

    def verify(self, symbolic_expr, domain, symbol_list):
        """Verify using dReal.

        Main wrapper function of the verification class.
        """
        # Initialize solver instance
        solver = z3.Solver()
        # Write SMT file
        string = self.make_SMT2_string(symbolic_expr,
                                       domain, symbol_list)
        smt_parse = z3.parse_smt2_string(string)
        solver.add(smt_parse)
        flag = solver.check()

        result = {}
        if str(flag) == "sat":
            # Create a list of z3 variables
            # (required to extract the counter examples)
            model = solver.model()
            result['sat'] = False
            z3_symbols = [z3.Real(str(i)) for i in symbol_list]
            result['violation'] = tuple([model[i].numerator_as_long(
            )/model[i].denominator_as_long() for i in z3_symbols])
        elif str(flag) == "unsat":
            result['sat'] = True
            result['violation'] = None
        else:
            result['sat'] = "unknown"
            result['violation'] = None

        return result

    def symbol_to_z3_symbol(self, symbol_list):
        """Convert a sympy Real symbol to a Z3 real symbol."""
        z3_symbols = [z3.Real(str(i)) for i in symbol_list]
        return z3_symbols


def demo():
    """Demo how to use this module."""
    import sympy as sp
    """Demonstration how to use the module."""
    # Define path variables
    path = 'e:/docker_connect'

    # Declare symbolic variables
    varlist = x, y, z = sp.symbols("x y z")
    # Specify domain
    dom = sp.And(x >= -1, x <= 1, y >= -1, y <= 1, z >= -1, z <= 1)
    # Declare symbolic inequality
    F = x*y+1 >= sp.cos(x/3.)+x/2.

    # dReal
    options = {'path': path}
    verifier = Dreal(options)
    result = verifier.verify(F, dom, varlist)
    print(result)

    # Z3
    F = x*y+1 >= x/3.+x/2
    verifier = Z3()
    result = verifier.verify(F, dom, varlist)
    print(result)

    return 1


# demo()
