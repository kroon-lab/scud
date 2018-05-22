#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

from scitbx.array_family import flex
from iotbx import mtz

import phil as sigma2mtz_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.mtz.MTZ import mtz_correlation

def run(args=None, l= None):
    '''
    add Sigma columns to output from scud.pdb.supercell
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='sigma2mtz_log.txt')
    l.title("exp.sigma2mtz module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = sigma2mtz_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.sigma2mtz

###########################################################################
#                           Read mtz_files                                #
###########################################################################

    l.show_info('Processing mtz file: {}'.format(p.input.mtz_in))
    # Read mtz, create mtz_object
    mtz_object = mtz.object(p.input.mtz_in)
    # Extract miller arrays from mtz_object
    miller_arrays = mtz_object.as_miller_arrays()
    # create new miller array (IBRG) with added sigma's (sqrt(IBRG))
    ibrg_new = miller_arrays[0].customized_copy(data=miller_arrays[0].data(),
                                                sigmas=flex.sqrt(miller_arrays[0].data()))
    # Create new mtz dataset
    mtz_dataset = ibrg_new.as_mtz_dataset(column_root_label="IBRG")
    # create new miller array (ITOT) with added sigma's (sqrt(ITOT))
    itot_new = miller_arrays[1].customized_copy(data=miller_arrays[1].data(),
                                                sigmas=flex.sqrt(miller_arrays[1].data()))
    mtz_dataset.add_miller_array(itot_new,
                                 column_root_label="ITOT")
    # create new miller array (IDFF) with added sigma's (sqrt(IDFF))
    idff_new = miller_arrays[1].customized_copy(data=miller_arrays[2].data(),
                                                sigmas=flex.sqrt(miller_arrays[2].data()))
    mtz_dataset.add_miller_array(idff_new,
                                 column_root_label="IDFF")
    
    # Create miller array from ITOT with lower sigma's on non-Bragg
    # Works only for 222 right now

    n = p.params.n_sigma
    
    # Get itot array
    itot2 = miller_arrays[1]
    # Select indices and create mask for uneven values
    itot2_indices = np.array(itot2.indices())
    m = flex.bool(((itot2_indices % 2).sum(axis=1) != 0))
    
    # Create sqrt(itot) array
    new_sig = flex.sqrt(itot2.data())
    # Set uneven miller index sigma's to 1/100 * sigma
    new_sig.set_selected(m, new_sig.select(m) / 100)

    if p.params.IDFF_none == False:
        d = itot2.data()
    if p.params.IDFF_none == True:
        d = itot2.data()
        d.set_selected(m , False)


    # new miller array with new ITOT and ITOT:SIGMA
    itot2_new = miller_arrays[1].customized_copy(data=d,
                                                 sigmas=new_sig)
    # Add to mtz_dataset
    mtz_dataset.add_miller_array(itot2_new,
                                 column_root_label="ITOT2")
    # Write new MTZ to file
    mtz_dataset.mtz_object().write("{}_SIG.mtz".format(p.input.mtz_in[:-4]))

    
