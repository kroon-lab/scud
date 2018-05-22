#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os
import random

from scitbx.array_family import flex

import phil as simple_trim_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass

def run(args=None, l= None):
    '''
    Remove nr models from large PDB file
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='simple_trim_log.txt')
    l.title("exp.simple_trim module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = simple_trim_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.simple_trim

###########################################################################
#                      Loop over all lines                                #
###########################################################################

    # Running parameters
    w = True
    c = p.params.nr
    cnt = 0

    l.show_info('Skipping {} models'.format(c))
    
    # Loop over lines in input PDB
    with open(p.input.pdb_in,'r') as f:
        with open(p.input.pdb_in[:-4]+'_shrunk.pdb','w') as n:
            for line in f:
                # Write if w = True
                if w == True:
                    print >> n, line.rstrip()
                # Check if endmdl line and cnt
                if line[:6] == 'ENDMDL':
                    if (cnt % c) == 0:
                        w = True
                    else:
                        w = False
                    cnt += 1
