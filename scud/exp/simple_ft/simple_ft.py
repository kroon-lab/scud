#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import mmtbx.utils
from libtbx import easy_mp
from iotbx import mtz

from scitbx.array_family import flex

import phil as simple_ft_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass
import scud.pdb.supercell.supercell as supercell

def run(args=None, l= None):
    '''
    Check additivity of motion in reciprocal space, subtract mtz files
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='additivity_log.txt')
    l.title("exp.simple_ft module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = simple_ft_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.simple_ft

    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    mtz_ls = ['tmp_{}.mtz'.format(i) for i in range(len(pdb_object.hierarchy.models()))]
    print mtz_ls
    pdb_object.set_B_zero()

    if p.params.FFT == True:
        pdb_object.start_or_attach_hierarchy()
        for i,mdl in enumerate(pdb_object.hierarchy.models()):
            # convert hierarchy to Xray structure
            tmp_model = pdb_object.detached_model_copy(n=i)
            pdb_object.start_or_attach_hierarchy(detached_model=tmp_model)
            xray_structure = pdb_object.new_hierarchy.extract_xray_structure(
                crystal_symmetry = pdb_object.symmetry)
            # Get default fmodel parameters
            params = mmtbx.command_line.fmodel.\
                     fmodel_from_xray_structure_master_params.extract()
            # Input high resolution
            params.high_resolution = 2.0
            mtz_name = 'tmp_{}.mtz'.format(i)
            mtz_ls.append(mtz_name)
            # Convert structure to structurefactors
            mmtbx.utils.fmodel_from_xray_structure(
                xray_structure = xray_structure,
                f_obs          = None,
                add_sigmas     = False,
                params         = params,
                twin_law       = None,
                twin_fraction  = None,
                out            = sys.stdout).write_to_file(file_name = mtz_name,
                                                           obs_type=complex)

    for c,mtz_name in enumerate(mtz_ls):
        print c, mtz_name
        if c == 0:
            template_object = mtz.object(mtz_name)
            template_miller_arrays = template_object.as_miller_arrays()
            f = template_miller_arrays[0]
            i = f.as_intensity_array().data()
        else:
            mtz_object = mtz.object(mtz_name)
            miller_arrays = mtz_object.as_miller_arrays()
            f_temp = miller_arrays[0]
            i_temp = f_temp.as_intensity_array().data()
            f += f_temp
            i += i_temp
    # Average f and i and return
    f_average = f / len(args)
    i_average = i / len(args)
    i_brg = f_average.as_intensity_array().set_observation_type_xray_intensity()
    # Calculate i_diff by I_tot - I_brg
    i_diff = i_average - i_brg.data()
    # Prep miller arrays for output using i_brg as template
    i = i_brg.customized_copy(data=i_average)
    i_diff = i_brg.customized_copy(data=i_diff)
    # Create mtz dataset, add miller_arrays with proper column labels.
    mtz_dataset = i_brg.as_mtz_dataset(column_root_label="IBRG",column_types='J')
    mtz_dataset.add_miller_array(i,column_root_label="ITOT")
    mtz_dataset.add_miller_array(i_diff,column_root_label="IDFF")
    # Write to file
    mtz_dataset.mtz_object().write("result.mtz")
    
    return l
    
