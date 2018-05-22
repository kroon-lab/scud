#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import iotbx.pdb
from scitbx.array_family import flex

import phil as cypa_helix_phil
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
        l = Log(log_name='cypa_helix_log.txt')
    l.title("exp.cypa_helix module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = cypa_helix_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.cypa_helix

###########################################################################
#                           Read PDB file                                 #
###########################################################################

    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")

    # HELIX residues: 30-42
    res_list = range(30,43)
    t_list = np.arange(-.50,0.52,0.01)
    final_hierarchy = iotbx.pdb.hierarchy.root()
    for t in t_list:
        new_model = pdb_object.hierarchy.models()[0].detached_copy()
        for resi in res_list:
            for at in new_model.chains()[0].residues()[resi].atoms():
                at.xyz = tuple(np.array(at.xyz)+np.array([0,0,t]))
        final_hierarchy.append_model(new_model)

    final_hierarchy.write_pdb_file(crystal_symmetry = pdb_object.symmetry,
                                   file_name = p.input.pdb_in[:-4]+'_helix.pdb')
        
