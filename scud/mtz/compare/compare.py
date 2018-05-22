#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as compare_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.mtz.MTZ import mtz_correlation

def run(args=None):
    '''
    Calculates R-value and CC values between 2 miller arrays 
    supplied in .mtz format
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='compare_log.txt')
    l.title("mtz.compare module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = compare_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.compare

###########################################################################
#                         Reading .mtz files                              #
###########################################################################

    # Check if mtz files are present
    if p.input.mtz_1 != None and p.input.mtz_2 != None:
        l.process_message('Reading mtz files...')
        mtz_1 = MTZClass(p.input.mtz_1,p.params.array_1)
        l.show_info('mtz file 1 {} read'.format(p.input.mtz_1))
        mtz_2 = MTZClass(p.input.mtz_2,p.params.array_2)
        l.show_info('mtz file 2 {} read'.format(p.input.mtz_2))
    else:
        l.warning('Not enough mtz files inputted, exiting')
        exit()

    l.process_message('Calculating R-value and CC...')
    mtz_correlation(mtz_1 = mtz_1,
                    mtz_2= mtz_2,
                    super_1 = p.params.supercell_size_1,
                    super_2 = p.params.supercell_size_2,
                    l = l)

###########################################################################
#                            Close Log File                               #
###########################################################################

    l.process_message('Finished, exiting...')
    l.close_log()
