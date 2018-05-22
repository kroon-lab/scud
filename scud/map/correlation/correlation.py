#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import uctbx
import iotbx.ccp4_map
from cctbx import sgtbx

import phil as correlation_phil
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
    l = Log(log_name='correlation_log.txt')
    l.title("map.correlation module")

###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = correlation_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.correlation

###########################################################################
#                              Read Map                                   #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_1))
    map_1 = MAPClass(p.input.map_1)
    l.process_message('Opening map {}'.format(p.input.map_2))
    map_2 = MAPClass(p.input.map_2)

###########################################################################
#                            Print Info                                   #
###########################################################################

    l.warning('map 1 info: ')
    map_1.map_info(l=l)
    l.warning('map 2 info: ')
    map_2.map_info(l=l)

    # Map info
    grid_1 = map_1.maps.unit_cell_grid
    grid_2 = map_2.maps.unit_cell_grid
    uc_1 = map_1.uc
    uc_2 = map_2.uc
    grd = flex.grid(grid_1)
    
    dat_1 = map_1.maps.data.deep_copy()
    dat_2 = map_2.maps.data.deep_copy()

    print dat_1.as_1d().size()
    
    # Find filter
    if p.params.use_all:
        values = ( dat_1  > -10e14 )
    if p.params.filter_map_1:
        values = ( dat_1.as_double().as_1d() != -1000 )
        t2 =  ( dat_2.as_double().as_1d() != -1000 )
        #assert sum(values) != values.size()

    print dat_2.as_double().as_1d().select(t2).size()
        
    # Filter both (1d) lists
    d1_filtered = dat_1.as_double().as_1d().select(values)
    d2_filtered = dat_2.as_double().as_1d().select(values)*3.5 + 14.1
 
    # print info
    print d1_filtered.size()

    # Calc cc
    mc = flex.linear_correlation(x = d1_filtered,
                                 y = d2_filtered).coefficient()

    l.show_info('Map Map correlation: {}'.format(mc))

    # Treat both maps the same and write to file (except m2 is filtered)
    other_val = ( dat_1.as_double().as_1d() == -1000 ) # Where is m1 -1000
    dat_2_f = dat_2.as_double().as_1d().set_selected(other_val,-1000) # set m2 -1000 based on other val
    dat_2_f.reshape(grd)
    
    d1 = dat_1.as_double().as_1d()
    d1.reshape(grd)

    # Write maps to file
    print 'm1'
    iotbx.ccp4_map.write_ccp4_map(        
            file_name= 'm1_ref.map',
            unit_cell = uctbx.unit_cell(uc_1),
            space_group= sgtbx.space_group_info("P1").group(),
            map_data = d1,
            labels=flex.std_string(["iotbx.ccp4_map.tst"]))

    print 'm2'
    iotbx.ccp4_map.write_ccp4_map(        
            file_name= 'm2_ref.map',
            unit_cell = uctbx.unit_cell(uc_1),
            space_group= sgtbx.space_group_info("P1").group(),
            map_data = dat_2_f,
            labels=flex.std_string(["iotbx.ccp4_map.tst"]))
