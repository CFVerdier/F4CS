# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 17:36:11 2021

@author: Cees F. Verdier
"""

from setuptools import setup

setup(
    name='f4cs',
    version='0.1',
    author='Cees Ferdinand Verdier',
    license='LICENSE.md',
    description=('A tool for (Formal) Correct-by-Construction Closed-form '
                 'Controller Synthesis'),
    long_description=open('README.md').read(),
    install_requires=[
        'numpy',
        'sympy',
        'cma',
        'z3-solver',
        'pyibex',
        'docker ; platform_system=="Windows"'
    ],
    python_requires=">=3.7"
)
