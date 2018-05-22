#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as slice_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.general.Plot import Plot
from scud.general.method_library import input_mtz
from scud.general.method_library import datetime_string

def run(args=None):
    '''
    Extract l=0 slice from mtz file, export as new .mtz, .npy and plot
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='slice_mtz_log.txt')
    l.title("exp.check_2d_detail")
    
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
        mtz_obj = MTZClass(p.input.mtz_in)

###########################################################################
#                        Slicing .mtz file                                #
###########################################################################

        # Select which array to slice
        l.process_message('Slicing mtz file and writing to {}...'.format(p.output.mtz_out))
        ar = p.params.array
        if ar == 'Ibrg':
            ar_num = 0
            l.warning('Array name: Ibrg')
        if ar == 'Itot':
            ar_num = 1
            l.warning('Array name: Itot')
        if ar == 'Idff':
            ar_num = 2
            l.warning('Array name: Idff')
        slice_2D = mtz_obj.slice_mtz(mtz_out_name=p.output.mtz_out,array_num=ar_num)
        l.process_message('Saving 2D slice as .npy file: {}...'.format(p.output.npy_out))
        np.save(p.output.npy_out, slice_2D)

###########################################################################
#                          Reading .npy Slice                             #
###########################################################################

    if p.input.npy_in != None:
        l.process_message('Loading 2D array from npy file...')
        slice_2D = np.load(p.input.npy_in)

###########################################################################
#                          Cut 2D map!                                    #
###########################################################################

    sp = np.shape(slice_2D)[0]+1
    delta = int(round(sp/10))
    delta_y = int(round(delta/2))
    center = int((sp / 2) - 1)

    print sp, delta, delta_y, center

    slice_2D = slice_2D[-delta:,center-delta_y:center+delta_y]

    print np.shape(slice_2D)

###########################################################################
#                             Plotting Slice                              #
###########################################################################

    l.process_message('Plotting 2D slice and saving to: {}...'.format(p.output.plt_out))
    l.show_info('Minimum of slice: {:.2f}'.format(np.nanmin(slice_2D)))
    l.show_info('Maximum of slice: {:.2f}'.format(np.nanmax(slice_2D)))
    if p.params.vmin == None and p.params.vmax == None:
        l.warning('Automatic determination of minimum and maximum values')
    else:
        l.show_info('Minimum value for plotting: {:.2f}'.format(p.params.vmin))
        l.show_info('Maximum value for plotting: {:.2f}'.format(p.params.vmax))        
    pt = Plot()
    pt.contour2D(slice_2D,plotfile_name=p.output.plt_out,vmin=p.params.vmin,vmax=p.params.vmax)
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
