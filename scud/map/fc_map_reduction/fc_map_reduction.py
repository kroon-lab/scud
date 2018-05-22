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

# MTZ stuff
from iotbx import mtz
from cctbx import miller
from cctbx import crystal

import phil as fc_map_reduction_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass

#from scud.mtz.MTZ import MTZClass

def run(args=None):
    '''
    - read mtz file, select F+phi and convert to ED map (.map)
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='fc_map_reduction_log.txt')
    l.title("map.fc_map_reduction module")

###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = fc_map_reduction_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.fc_map_reduction

###########################################################################
#                              Read MTZ                                   #
###########################################################################

    # Should be moved to mtz/MTZ.py :
    l.process_message('Reading MTZ and extracting columns...')
    # Create mtz object
    mtz_object = mtz.object(p.input.mtz_in)
    # Hardcoded labels!!!
    label = 'F-model_xray'
    phi_lab = 'PHIF-model_xray'

    # Create custom miller set:
    # Crystal information from 'label':
    crystal = mtz_object.get_column(label).mtz_crystal().crystal_symmetry()
    # Extract amplitudes and phases
    fmodel_complex = mtz_object.extract_complex(column_label_ampl=label,  # amplitude label
                                                column_label_phi=phi_lab) # phase label
    # Actually create miller set
    miller_set = miller.set(crystal_symmetry=crystal,         # Crystal symmetry
                            indices=fmodel_complex.indices,   # Miller indices
                            anomalous_flag=False)             # No anomalous data
    # Create complex miller array
    miller_array = miller_set.array(fmodel_complex.data)

###########################################################################
#                           SF+PH -> ARRAY                                #
###########################################################################

    l.process_message('Creating electron density map and testing resolution factors...')

    # Create array form ed map, and check if dimensions are even
    # Start values:
    is_even = False
    start_rf = 0.33
    # loop:
    while is_even == False:
        l.show_info('Resolution factor: {}'.format(start_rf))
        # FFT amplitudes + phases to electron density (ed) map
        ed_map = miller_array.fft_map(resolution_factor=start_rf)
        # Create array from fft
        ed_array = ed_map.real_map_unpadded()
        # dimensions:
        dims = np.array(ed_array.all())
        # Check if dimensions are  even
        if np.sum(dims % 2) != 0:
            start_rf -= 0.01
        else:
            is_even = True
        
    # Write whole ed map to ccp4 type map
    ed_map.as_ccp4_map(file_name = 'supercell_ed.map')

    l.show_info('Supercell Array size: {}'.format(ed_array.all()))
    l.process_message('Extracting unit cells from map...')

###########################################################################
#                           Array Slicing                                 #
###########################################################################

    # Slice array in 8 sub arrays and write to .map files
    l.process_message('Slicing array...')

    # Slicing parameters
    limits = ((dims / 2)-1).astype(int)
    print limits
    m_list = map_slicer(ed_array, limits)


    # New unit cell (same for all)
    # super cell unit cell:
    sc_uc = crystal.unit_cell().parameters()
    new_uc = uctbx.unit_cell((sc_uc[0]/2,
                              sc_uc[1]/2,
                              sc_uc[2]/2,
                              sc_uc[3],
                              sc_uc[4],
                              sc_uc[5]))

    # Write maps:
    for i,m in enumerate(m_list):
        print m.all()
        nm = 'm{}.map'.format(i)
        iotbx.ccp4_map.write_ccp4_map(        
            file_name= nm,
            unit_cell = new_uc,
            space_group= sgtbx.space_group_info("P1").group(),
            map_data = m,
            labels=flex.std_string(["iotbx.ccp4_map.tst"]))

###########################################################################
#                              P1 -> ASU
#
###########################################################################
    

    return l

def map_slicer(ar, l):
    
    m1 = ar[0:l[0],
            0:l[1],
            0:l[2]]
    m2 = ar[l[0]+1:-1,
            0:l[1],
            0:l[2]]
    m3 = ar[0:l[0],
            l[1]+1:-1,
            0:l[2]]
    m4 = ar[0:l[0],
            0:l[1],
            l[2]+1:-1]
    m5 = ar[l[0]+1:-1,
            l[1]+1:-1,
            0:l[2]]
    m6 = ar[l[0]+1:-1,
            0:l[1],
            l[2]+1:-1]
    m7 = ar[0:l[0],
            l[1]+1:-1,
            l[2]+1:-1]
    m8 = ar[l[0]+1:-1,
            l[1]+1:-1,
            l[2]+1:-1]
    return [m1,m2,m3,m4,m5,m6,m7,m8]
    
