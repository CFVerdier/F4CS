# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 12:59:31 2021

@author: ceesv
"""

import subprocess


def call_dReal(specification):
    """Call dReal.

    OS specific
    """
    if specification.t_max is None:
        outputdReal = subprocess.check_output(
            ['powershell.exe',
             "docker run -v " + specification.path
             + ":/data --rm dreal/dreal4 dreal data/" + specification.file_name
             + " --model"],
            shell=True, stderr=subprocess.STDOUT).decode("utf-8")
    else:
        outputdReal = subprocess.check_output(
            ['powershell.exe',
             "docker run -v " + specification.path
             + ":/data --rm dreal/dreal4 dreal data/" + specification.file_name
             + " --model"],
            shell=True,
            stderr=subprocess.STDOUT,
            timeout=specification.t_max).decode("utf-8")

    return outputdReal
