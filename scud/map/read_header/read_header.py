#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as read_header_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass

def run(args=None):
    '''
    Read a ccp4 type map and output information
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='read_header_log.txt')
    l.title("map.read_header module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = read_header_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.read_header

###########################################################################
#                              Read Map                                   #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_in))
    map_obj = MAPClass(p.input.map_in)

###########################################################################
#                            Print Info                                   #
###########################################################################

    map_obj.map_info(l=l)

###########################################################################
#                         Shift Origin                                    #
###########################################################################

    if p.params.shift_origin == True:
        new_map = map_obj.shift_origin(l=l)
    new_obj = MAPClass(new_map)
    new_obj.map_info(l=l)
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
