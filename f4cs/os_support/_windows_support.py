# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 12:59:31 2021

@author: ceesv
"""

import docker


def call_dReal(verifier):
    """Call dReal.

    OS specific
    """
    path = verifier.path
    file = verifier.file_name
    client = docker.from_env()
    volume_dict = {path: {'bind': '/data', 'mode': 'ro'}}
    container = client.containers.run("dreal/dreal4",
                                      " dreal /data/" + file + " --model",
                                      volumes=volume_dict,
                                      detach=True)
    if verifier.t_max is not None:
        container.stop(timeout=verifier.t_max)
    outputdReal = container.logs().decode("utf-8")
    container.remove()
    if outputdReal == '':
        print("dReal time out: {}".format(verifier.t_max))
        outputdReal = 'time-out'

    return outputdReal
