#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as slice_map_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass
from scud.plt.Plot import Plot

def run(args=None):
    '''
    slice and plot hk0, h0l, and 0kl slices from map and plot
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='slice_map_log.txt')
    l.title("map.slice_map module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = slice_map_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.slice_map

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
#                            Slice Map                                    #
###########################################################################

    # shape of the map
    shape = map_obj.maps.data.all()
    # Origin of map
    origin = map_obj.origin_shift
    # Data in 3D array
    data = map_obj.maps.data.as_numpy_array()
    l.show_info('Data min: {}'.format(np.min(data)))
    l.show_info('Data max: {}'.format(np.max(data)))
    # Slice hk0, h0l, 0kl planes
    data_1 = np.array(data[origin[0],
                           :,
                           :])
    data_2 = np.array(data[:,
                           origin[1],
                           :])
    data_3 = np.array(data[:,
                           :,
                           origin[2]])

    m_val = data_1[0,0]
    
    data_1[(data_1 == m_val)] = np.nan
    data_2[(data_2 == m_val)] = np.nan
    data_3[(data_3 == m_val)] = np.nan
    
###########################################################################
#                            Plot Slices                                  #
###########################################################################


    mp_ls = [data_1, data_2, data_3]
    prefix = p.input.map_in[:-4]
    nms = [prefix+'_h.eps',prefix+'_k.eps',prefix+'_l.eps']

    vmin = p.params.vmin
    vmax = p.params.vmax
    
    for i,dat in enumerate(mp_ls):
        pt = Plot()
        pt.contour2D(dat,
                     plotfile_name=nms[i],
                     vmin=vmin, vmax=vmax,
                     cnt_num=15)
    pt = Plot()
    pt.colorbar(nm=p.input.map_in[:-4]+'_hkl0_cmap.eps',
                vmin = vmin,
                vmax = vmax)

###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
