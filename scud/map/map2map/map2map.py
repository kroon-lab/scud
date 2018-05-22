#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os


from scitbx.array_family import flex

import phil as map2map_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass
from scud.map.MAP import MapMapClass
from scud.plt.Plot import Plot

def run(args=None,
        l = None):
    '''
    * Read 2 ccp4 type maps. 
    * Cut the larger one to be similar to the smaller map
    * Shift both origins to 0,0,0    
    * Write both 'new maps'
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='map2map_log.txt')
    l.title("map.map2map module")
    # Plot class:
    pl = Plot()
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.warning('If one of the maps is a data map, input it as map 1 and set map_1_exp=True')
    l.process_message('Processing input...\n')
    # use phil to process input
    working_params = map2map_phil.phil_parse(args=args,log=l)  
    p = working_params.map2map

###########################################################################
#                              Read Maps                                  #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_1))
    map_1 = MAPClass(p.input.map_1)
    l.process_message('Opening map {}'.format(p.input.map_2))
    map_2 = MAPClass(p.input.map_2)

    if p.params.map_1_exp: l.show_info('{} is expermental map'.format(p.input.map_1))
    if p.params.map_2_exp: l.show_info('{} is expermental map'.format(p.input.map_2))
        
###########################################################################
#                            Print Info                                   #
###########################################################################

    l.process_message('map 1 ({}):'.format(p.input.map_1))
    map_1.map_info(l=l)
    l.show_info('\n')
    l.process_message('map 2 ({}):'.format(p.input.map_2))
    map_2.map_info(l=l)

###########################################################################
#                        Do Map Map analysis                              #
###########################################################################

    ######!!!!!#######
    # TODO:
    # Make Mask on non-Brag intensities, calculate CC on Bragg
    # Radial intensity plots on IDFF

    l.show_info('\n')
    # Load MapMap class (contains map map analysis tools)
    mm = MapMapClass(map_1 = map_1,
                     map_1_dat = p.params.map_1_exp,
                     map_2 = map_2,
                     map_2_dat = p.params.map_2_exp,
                     p = p,
                     l = l)
    # Resize maps so both arrays are same size
    mm.resize_maps(l = l,
                   p = p)
    # Create mask based on padding values
    mm.create_mask(p = p,
                   l = l)
    # Calculate Cross Correlation between maps
    mm.calc_CC(p = p,
               l = l)
    # Scale maps to later calc difference map
    mm.scale_maps(p = p,
                  l = l)
    # Calculate radial profile of both maps and write plot
    mm.radial(p = p,
              l = l)
    # Calculate difference map(s)
    mm.diff_map(p = p,
                l = l)
    # Project difference map
    mm.project_diff(p = p,
                    l = l)

###########################################################################
#                             Close Log file                              #
###########################################################################

    return l
