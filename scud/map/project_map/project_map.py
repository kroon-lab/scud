#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as project_map_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass

def run(args=None):
    '''
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='project_map_log.txt')
    l.title("map.project_map module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    # use phil to process input
    working_params = project_map_phil.phil_parse(args=args,log=l)  
    p = working_params.project_map

###########################################################################
#                              Read Map                                   #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_in))
    map_in = MAPClass(p.input.map_in)
    map_in.map_info(l = l)
    # Project map
    map_in.projection(l = l,
                      outname = p.output.plot_out,
                      mask_value = p.params.mask_value)
    
###########################################################################
#                             Close Log file                              #
###########################################################################

    l.close_log()

