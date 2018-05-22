#!/usr/bin/env cctbx.python

import numpy as np
import sys,os

from scud.pdb.PDB import PDBClass
import phil as trim_ensemble_phil
from libtbx.phil import parse
from scud.general.log_writer import Log

def run(args=None, l=None):
    """
    Generate an ensemble describing types of rigid body motion,
    input can be single structure or ensemble describng 'internal' motion
    """
###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    if l == None:
        l = Log(log_name='trim_ensemble_log.txt')
    l.title("pdb.trim_ensemble module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    working_params = trim_ensemble_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.trim_ensemble
    
###########################################################################
#                           Trim ensemble                                 #
###########################################################################

    l.process_message('Reading and analyzing PDB file / Structure...')
    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")

    l.show_info("Number of models in ensemble: {}".format(pdb_object.model_num))
    l.show_info("Will be trimmed to {} models".format(p.params.ens_size))
    l = pdb_object.trim_ensemble(target = p.params.ens_size,
                                 out_name = p.output.pdb_out,
                                 skip_method = True,
                                 l = l)
    
    return l
