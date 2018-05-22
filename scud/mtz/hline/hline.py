#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os
import json

import phil as hline_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass

def run(args=None):
    '''
    Extract k=0 l=0 slice from mtz file, export as new .mtz, .npy and plot
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='hline_log.txt')
    l.title("mtz.hline module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = hline_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.hline

###########################################################################
#                         Reading .mtz file                               #
###########################################################################

    if p.input.mtz_in != None:
        l.process_message('Reading mtz file...')
        mtz_obj = MTZClass(p.input.mtz_in,p.params.array)

###########################################################################
#                        Slicing .mtz file                                #
###########################################################################

    l.process_message('Slicing mtz file')
    data = mtz_obj.hline_mtz(supercell_size=p.params.supercell_size)
    
    # normalize:
    h = [i[0] for i in data]
    d = np.array([i[1] for i in data])
    d_norm = d / np.amax(d)
    # output
    a = {}
    a['xlabel'] = 'h'
    a['ylabel'] = 'norm(I)'
    a['x'] = list(h)
    a['y'] = list(d_norm)
    a['label'] = str(p.params.supercell_size)
    a['ls'] = '-'
    json.dumps(a,sort_keys=True)
    with open('outfile.json','w') as outfile:
        json.dump(a,outfile)

    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
