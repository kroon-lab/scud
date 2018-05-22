#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import iotbx.ccp4_map
from cctbx.array_family import flex
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import crystal
from cctbx import miller
from iotbx import mtz

import phil as map2mtz_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.map.MAP import MAPClass

def run(args=None, l=None):
    '''
    Convert structure factor file (mtz) to ccp4 type map.
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='map2mtz_log.txt')
    l.title("mtz.map2mtz module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = map2mtz_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.map2mtz

###########################################################################
#                         Reading .mtz files                              #
###########################################################################

    map_in = MAPClass(fname=p.input.map_in)
    map_in.map_info(l = l)
    l.process_message('Map file read...')

###########################################################################
#                    Initializing final map                               #
###########################################################################

    # Data array as numpy array
    dat = map_in.maps.data.as_numpy_array()
    # Find minimum value (not unmeasured)
    if p.params.make_positive == True:
        f = np.sort((dat.flatten()))
        begin = f[0]
        for i in f:
            if i != begin:
                minimum = abs(i)
                break
        
        l.show_info('Minimum = {}'.format(minimum))
    # Origin shift (difference between index and hkl value)
    shift = np.array(map_in.origin_shift).astype(int) - 1
    # Shape of numpy array / data
    sp = np.shape(dat)
    # List to save hkl and intensity
    hkl_I = []
    # Loop over List (realy slow)
    l.process_message('(Re-)indexing Map...')
    minimum = 0.0
    for i in range(sp[0]):
        for j in range(sp[1]):
            # Start at l = 0 (centrosymmetry)
            for k in range(shift[2],sp[2]):
                if dat[i,j,k] != -1000.0:
                    ind = np.array([i,j,k])
                    hkl_I.append([tuple(ind-shift),dat[i,j,k]+minimum])

###########################################################################
#                           Create MTZ file                               #
###########################################################################

    l.process_message('Creating miller object...')
    # Create a flex miller list
    mi = flex.miller_index(zip(*hkl_I)[0])
    # Create symmetry object
    vxs = map_in.voxel_size
    angl = map_in.uc[3:7]
    symmetry = crystal.symmetry(unit_cell = (np.append(np.array(vxs),np.array(angl))),
                                space_group_symbol="P1")
    # Create miller set
    ms = miller.set(
        crystal_symmetry = symmetry,
        anomalous_flag = False,
        indices = mi)

    l.process_message('Converting data to flex...')
    # Convert data Intensities to flex double
    new_data = flex.double(zip(*hkl_I)[1])
    # Create miller array
    ma = miller.array(ms,
                      data=new_data)
    l.process_message('Creating mtz and writing to file...')
    # Create mtz_data set and write
    mtz_dataset = ma.as_mtz_dataset(column_root_label='IDFF',column_types='J')
    mtz_dataset.mtz_object().write(p.input.map_in[:-4]+'_new.mtz')

    return l
