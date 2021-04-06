# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 12:59:31 2021

@author: ceesv
"""

import subprocess
import os


def call_dReal(verifier):
    """Call dReal.

    OS specific.
    """
    dreal_path = verifier.dReal_path
    base_path = verifier.path
    file_name = verifier.file_name

    if verifier.t_max is None:
        outputdReal = subprocess.check_output(
            [dreal_path, '--model',
             os.path.join(base_path, file_name)]).decode("utf-8")
    else:
        outputdReal = subprocess.check_output(
            [dreal_path, '--model', os.path.join(base_path, file_name)],
            timeout=verifier.t_max).decode("utf-8")

    return outputdReal
