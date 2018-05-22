from __future__ import division
from cctbx import crystal
from cctbx import uctbx
from cctbx import miller
import iotbx.pdb
from iotbx import mtz
import mmtbx.utils
from libtbx import easy_mp

import numpy as np
from random import randint
import os,sys

import phil as correlated_supercell_phil
from scud.pdb.PDB import PDBClass
from scud.general.log_writer import Log

def create_correlated_sc(l=None):

    l.process_message('Creating supercell with internal correlations in the displacement...')
    # pseudo code:

    # calculate center of mass of asymmetric unit
    # Generate list with all symmetry operations to fill a supercell
    # Append to each symmetry operation the coordinates of the center of mass after applying symmetry operation
    # Randomize order of list (So every supercell will be unique)
    # For the fist element of the list:
        # Place molecule which has a random translation vector
        # Check which other positions are within a certain distance (correlation lenght)
        # For each position calculate change it is correlated e^-Gamma
        # If it is correlated apply similar translation
        # If not apply random translation (also slightly modulated with respect to neighbour)
        # Choos random new molecule and repeat !?!?!?!?!?

    

def run(args=None, l=None):
    '''
    Create supercells containing slightly translated with respect to ideal 
    position. A correlation in the displacement is implemented.
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='correlated_supercell_log.txt')
    l.title("scud.pdb.correlated_supercell module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = correlated_supercell_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.correlated_supercell
    ss_num=p.input.supercell_num

    return l

