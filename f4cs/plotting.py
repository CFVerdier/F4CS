# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 21:14:27 2021

@author: Cees F. Verdier
"""
import numpy as np
import sympy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def plotting(variables, system, I_list, t_max=10, number_trajectories=5):
    """Generic plotting of closed-loop system trajectories."""

    n = len(variables)
    # Run multiple simulations
    tvec = np.linspace(0, t_max, 500)
    y = [tvec]*number_trajectories
    # Make the symbolic equations of motion a numerical function.
    t = sp.symbols('t')
    f_fun = sp.lambdify([variables, t], sp.flatten(system))

    # Bounds for the initial condition
    bounds = np.array(list(map(interval_width, I_list)))

    for i in range(number_trajectories):
        # random initial condition

        y0 = bounds * np.random.randint(-1, 2, size=n)
        # Random disturbances
        y[i] = odeint(f_fun, y0, tvec)

    # Plotting
    fig0, axes = plt.subplots(n, 1, figsize=(8, 7), sharex=True)
    for i, ax in enumerate(axes):
        ax.plot(tvec, np.array([x[:, i] for x in y]).T)
        ax.set_ylabel(f'x{i}')
    axes[n-1].set_xlabel('time [s]')
    plt.show()


def interval_width(interval):
    """Returns half the width of an interval, given by a list"""
    return (interval[1]-interval[0])/2+(interval[1]+interval[0])/2
