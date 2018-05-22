#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as slice_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.plt.Plot import Plot
from scud.general.method_library import input_mtz
from scud.general.method_library import datetime_string

def run(args=None, l=None):
    '''
    Extract l=0 slice from mtz file, export as new .mtz, .npy and plot
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='slice_mtz_log.txt')
    l.title("mtz.slice module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = slice_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.slice_mtz

###########################################################################
#                         Reading .mtz file                               #
###########################################################################

    if p.input.mtz_in != None:
        l.process_message('Reading mtz file...')
        mtz_obj = MTZClass(p.input.mtz_in,p.params.array)

###########################################################################
#                             Plotting Slice                              #
###########################################################################

    l.process_message('Plotting 2D slice and saving to: {}...'.format(p.output.plt_out))
    prefix = p.input.mtz_in[:-4]
    nms = [prefix+'_h.eps',prefix+'_k.eps',prefix+'_l.eps']

    vmin = p.params.vmin
    vmax = p.params.vmax

    # Slice data, loop over indices
    for i,ind in enumerate(['h','k','l']):
        slice_2D = mtz_obj.slice_mtz(slice_index=ind)
        if p.params.square:
            slice_2D = slice_2D**2
        min_val = np.nanmin(slice_2D)
        max_val = np.nanmax(slice_2D)
        l.show_info('Minimum of slice: {:.2f}'.format(min_val))
        l.show_info('Maximum of slice: {:.2f}'.format(max_val))
        if p.params.write_slices == True:
            n = nms[i][:-4]+'.npy'
            np.save(n,slice_2D)
        pt = Plot()
        # Plot
        pt.contour2D(slice_2D,
                     plotfile_name=nms[i],
                     vmin=vmin,vmax=vmax)

    pt = Plot()
    pt.colorbar(nm=p.input.mtz_in[:-4]+'_hkl0_cmap.eps',
                vmin = vmin,
                vmax = vmax)
    return l
