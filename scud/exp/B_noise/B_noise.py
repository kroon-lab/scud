#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os
import random

from scitbx.array_family import flex

import phil as B_noise_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass

def run(args=None, l= None):
    '''
    Check additivity of motion in reciprocal space, subtract mtz files
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='B_noise_log.txt')
    l.title("exp.B_noise module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = B_noise_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.B_noise

    # Read PDB
    pdb_f = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    # Loop over all atoms
    for atom in pdb_f.hierarchy.models()[0].atoms():
        atom.b = random.uniform(0.0,0.1)
        # Write PDB to file
    pdb_f.write_hierarchy_to_pdb(output_hierarchy=pdb_f.hierarchy,
                                 out_name='{}_B_noise.pdb'.format(p.input.pdb_in[:-4]))
