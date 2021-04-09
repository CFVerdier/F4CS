# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 17:36:11 2021

@author: ceesv
"""

from setuptools import setup

setup(
    name='f4cs',
    author='Cees Ferdinand Verdier',
    license='LICENSE.md',
    description='A tool for (Formal) Correct-by-Construction Closed-form Controller Synthesis',
    long_description=open('README.md').read(),
    install_requires=[
        'numpy',
        'sympy',
        'cma',
        'os',
        'platform',
        'collections',
        'subprocess',
        'z3-solver',
        'pyibex'
    ],
    python_requires=">=3.8"
)
