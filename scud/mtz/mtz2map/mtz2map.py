#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import iotbx.ccp4_map
from scitbx.array_family import flex
from cctbx import uctbx
from cctbx import sgtbx

import phil as mtz2map_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass

def run(args=None, l=None):
    '''
    Convert structure factor file (mtz) to ccpp4 type map.
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='mtz2map_log.txt')
    l.title("mtz.mtz2map module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = mtz2map_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.mtz2map

###########################################################################
#                         Reading .mtz files                              #
###########################################################################

    mtz = MTZClass(p.input.mtz_in,p.params.array)
    l.process_message('Input MTZ file read and column extracted...')
    if p.params.centro == False:
        if p.params.array2 != None:
            mtz2 = MTZClass(p.input.mtz_in,p.params.array2)
        else:
            print 'Give name of second array'
            quit()
###########################################################################
#                    Initializing final map                               #
###########################################################################

    # Calculate map shape
    min_max = mtz.miller_array.min_max_indices()
    l.show_info('min/max indices {}'.format(str(min_max)))
    h_range = (min_max[1][0]-min_max[0][0])+1
    k_range = (min_max[1][1]-min_max[0][1])+1
    l_range = (min_max[1][2]-min_max[0][2])+1
    # Calculate offset
    origin = tuple([min_max[0][0]*-1,
                    min_max[0][1]*-1,
                    min_max[1][2]])
    l.show_info('Origin/correction: {}'.format(str(origin)))

###########################################################################
#                             Fill the Map                                #
###########################################################################

    l.process_message('Filling final array...')
    # Final shape of the array and creation of empty array
    final_shape = (h_range,k_range, ((l_range-1)*2)+1)    
    zero_array = np.zeros(final_shape)
    # Get miller array, indices, centro-indices and their ravelled counterparts
    ma = mtz.miller_array
    idxs = np.array(ma.set().indices())
    vals = np.array(ma.data())
    if p.params.centro == False:
        ma2 = mtz2.miller_array
        vals2 = np.array(ma2.data())
    print p.params.array
    if p.params.array2 != None:
        print p.params.array2
    g = flex.grid(final_shape)
    # Mapping 3D miller indices to 1D array indices (slow)
    pos_idx = map(g, idxs + origin)
    neg_idx = map(g, (idxs*-1) + origin)
    # Fill the array (as 1D)
    if p.params.centro == True:
       zero_array.put(pos_idx, vals)
       zero_array.put(neg_idx, vals)
    else:
        zero_array.put(pos_idx,vals)
        zero_array.put(neg_idx,vals2)
    l.process_message('Filled array...')

###########################################################################
#                            Writing Map                                  #
###########################################################################

    l.process_message('Writing map to file...')
    # Final values for map creation (angles and length)
    angles = list(mtz.uc.parameters())[3:6]
    abc = mtz.uc.fractionalize(final_shape)
    l.show_info('a*b*c* map: {}'.format(str(abc)))
    # convert to flex format
    zero_array = flex.double(zero_array)
    # write map
    iotbx.ccp4_map.write_ccp4_map(        
        file_name = p.output.map_out,
        unit_cell = uctbx.unit_cell((abc[0],abc[1],abc[2],angles[0],angles[1],angles[2])),
        space_group = sgtbx.space_group_info("P1").group(),
        map_data = zero_array,
        labels = flex.std_string(["iotbx.ccp4_map.tst"]))
        
    l.process_message('Finished, exiting...')

    return l
