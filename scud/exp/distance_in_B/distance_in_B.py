#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os
import random

from scitbx.array_family import flex

import phil as distance_in_B_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass

def run(args=None, l= None):
    '''
    Compare ref PDB with input PDB
    calculate distance between center of masses
    Place in B-factor column of new PDB (copy of ref)
    Write PDB
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='distance_in_B_log.txt')
    l.title("exp.distance_in_B module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = distance_in_B_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.distance_in_B

    # reference (perfect) PDB
    ref_pdb = PDBClass(fname=p.input.ref_pdb, selection="not (resname HOH) and not (resname CL) and not (resname NA)")

    # test PDB
    t_pdb = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")

    tot_dists = []
    for ref_chain in ref_pdb.hierarchy.models()[0].chains():
        ref_com = ref_chain.atoms().extract_xyz().mean()
        d_list = []
        for t_chain in t_pdb.hierarchy.models()[0].chains():
            t_com = t_chain.atoms().extract_xyz().mean()
            d_list.append(np.linalg.norm(np.array(ref_com)-np.array(t_com)))
        new_b = np.min(d_list)
        for atom in ref_chain.atoms():
            atom.b = new_b
        tot_dists.append(np.min(d_list))

    ref_pdb.write_hierarchy_to_pdb(output_hierarchy=ref_pdb.hierarchy,
                                   out_name='dist_sc.pdb')

    import matplotlib.pyplot as plt
    plt.hist(tot_dists,
             bins = 10)
    plt.show()
