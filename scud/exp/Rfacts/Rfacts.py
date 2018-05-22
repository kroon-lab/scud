#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import mmtbx.utils
from libtbx import easy_mp
from iotbx import mtz

from scitbx.array_family import flex

import phil as Rfacts_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass
import scud.pdb.supercell.supercell as supercell
from scud.mtz.MTZ import MTZClass
from scud.mtz.MTZ import mtz_correlation

def run(args=None, l= None):
    '''
    First FFT structure
    Then calculate rvalues between output and other MTZ
    while splitting Itot and Idiff miller indices
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='Rfacts_log.txt')
    l.title("exp.Rfacts module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = Rfacts_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.Rfacts

###########################################################################
#                         FFT structure                                   #
###########################################################################

    # If input structure and data are to be compared FFT input structure first
    if p.input.er_mtz == None:
        pdb_1 = PDBClass(fname=p.input.pdb_1, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # create xray structure from model 1 in pdb_1
        xray_structure = pdb_1.hierarchy.extract_xray_structure(
            crystal_symmetry = pdb_1.symmetry)    
        # Get default fmodel parameters
        params = mmtbx.command_line.fmodel.\
                 fmodel_from_xray_structure_master_params.extract()
        # Input high resolution
        params.high_resolution = 2.0
        mtz_name = '{}.mtz'.format(p.input.pdb_1[:-4])
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
        l.process_message('Structure FFTed and written to file')

###########################################################################
#                   Read mtz's and calc R-factors                         #
###########################################################################

        # Read miller array 1
        mtz_1 = MTZClass('{}.mtz'.format(p.input.pdb_1[:-4]),p.params.array_1)
        # Read miller array 2
        mtz_2 = MTZClass(p.input.mtz_2,p.params.array_2)

    if p.input.er_mtz != None:
        array_1 = 'F-obs-filtered_xray'
        array_2 = 'F-model_xray'
        # Read er_mtz twice (2 columns needed)
        mtz_1 = MTZClass(p.input.er_mtz,array_1)
        mtz_2 = MTZClass(p.input.er_mtz,array_2)

    l.process_message('MTZ files read...')        
    l.process_message('Splitting data in IDFF and IBRG...')

    # Get miller array
    ma_1 = mtz_1.miller_array
    ma_2 = mtz_2.miller_array

    # Calculate scale factor between 2 columns
    sfactor = sf(ma_1 = ma_1,
                 ma_2 = ma_2,
                 l = l)
    # Mask Diffuse 
    ma_1_indices = np.array(ma_1.indices())
    m_brg = flex.bool(((ma_1_indices % 2).sum(axis=1) == 0))

    # Save R-values in this
    final_r = []
    # Calc R-Bragg
    nr = calculate_R(ma_1 = ma_1,
                     ma_2 = ma_2,
                     mask = m_brg,
                     name = 'Bragg',
                     sf = sfactor,
                     l = l)
    final_r.append(nr)
    
    # Mask Bragg
    m_dff = flex.bool(((ma_1_indices % 2).sum(axis=1) != 0))
    
    nr = calculate_R(ma_1 = ma_1,
                     ma_2 = ma_2,
                     mask = m_dff,
                     name = 'IDFF',
                     sf = sfactor,
                     l = l)
    final_r.append(nr)

    # Write to File
    if p.input.er_mtz != None:
        nm = '{}_ER_Rval.txt'.format(p.input.er_mtz[:-4])
    else:
        nm = '{}_input_Rval'.format(p.input.mtz_2[:-4])
    with open(nm,'w') as f:
        for i in final_r:
            print >> f, i
    
    return l

def sf(ma_1=None,
       ma_2=None,
       l = None):
    '''
    Calculate scale factor between all reflections (Br and Df)
    '''
    # Mask all (nothing)
    m = flex.bool(np.ones(ma_1.size()))
    sel_1 = ma_1.select(m)    
    select_ma_1 = ma_1.common_set(sel_1)
    select_ma_2 = ma_2.common_set(sel_1)
    # Convert intensities to amplitudes
    l.process_message('Converting to complex amplitude array')
    select_ma_1_amp = select_ma_1.as_amplitude_array().set_observation_type_xray_amplitude()
    select_ma_2_amp = select_ma_2.as_amplitude_array().set_observation_type_xray_amplitude()
    return select_ma_1_amp.scale_factor(select_ma_2_amp.randomize_phases())
    
def calculate_R(ma_1 = None,
                ma_2 = None,
                mask = None,
                name = None,
                sf = None,
                l = None):
    
    # Make sure both sets have the same reflections
    l.process_message('Selecting similar miller indices...')
    sel_1 = ma_1.select(mask)    
    select_ma_1 = ma_1.common_set(sel_1)
    select_ma_2 = ma_2.common_set(sel_1)
    # Convert intensities to amplitudes
    l.process_message('Converting to complex amplitude array')
    select_ma_1_amp = select_ma_1.as_amplitude_array().set_observation_type_xray_amplitude()
    select_ma_2_amp = select_ma_2.as_amplitude_array().set_observation_type_xray_amplitude()
    l.process_message('Calculating R value...')
    # Calculate R-value
    r_value = select_ma_1_amp.r1_factor(select_ma_2_amp,
                                        scale_factor=sf)
    l.show_info('{} R-value: {:.2f}'.format(name,r_value))
    return [name,r_value]
    
