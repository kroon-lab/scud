from __future__ import division
from cctbx import crystal
from cctbx import uctbx
from cctbx import miller
import iotbx.pdb
from iotbx import mtz
import mmtbx.utils
from libtbx import easy_mp
from cctbx.array_family import flex

import numpy as np
import os,sys
import random

import phil as sc_er_setup_phil
from scud.general.log_writer import Log

from scud.pdb.PDB import PDBClass
from scud.pdb.supercell import supercell
from scud.pdb.rigid_body import rigid_body

from scud.exp.EXP import EXPClass

def run(args=None, l=None):
    '''
    Create 4 PDB files, 2 P1, 2 SC, perfect or containing an offset
    Calculate (diffuse) scattering from all 4
    Run ensemble refinement for 8 combo's (4 per "size" SC or P1)
    
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='sc_er_setup_log.txt')
    l.title("er.sc_er_setup module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = sc_er_setup_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.sc_er_setup

    single_pdb_name = 'single_pdb.pdb'
    mtz_list, pdb_list = [],[]
    
###########################################################################
#                      Create simple dynamic model                        #
###########################################################################

    # Rigid Body motion simple test:
    if p.params.make_rb:
        l.process_message('RB module...')
        rigid_body_args = ['pdb_in={}'.format(p.input.rb_pdb_in),
                           'rb_type={}'.format(p.params.rb_type),
                           'pdb_out={}'.format(p.output.rb_pdb_out)]
        rigid_body.run(args = rigid_body_args,
                       l = l)
        os.system("mv trans_rb.pdb rb_pdb_out.pdb")
    
###########################################################################
#                          Prepare single PDB                             #
###########################################################################

    # Single PDB will be used for rb setups
    pdb_list.append(single_pdb_name)
    if p.params.prep_single:
        l.process_message('Prepping single_pdb')
        single_pdb = PDBClass(fname=p.input.single_pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Set B-factors 
        single_pdb.set_B_zero()
        # Set Occupancy
        single_pdb.set_occ()
        # Write PDB
        single_pdb.write_hierarchy_to_pdb(out_name=single_pdb_name,
                                          output_hierarchy=single_pdb.hierarchy)
   
###########################################################################
#                Calculate perfect scattering from single_pdb             #
###########################################################################

    if p.params.single_sc:
        l.process_message('single_pdb supercell operation...')
        # use EXP class to calculate diffuse scattering
        ex = EXPClass()
        ex.calc_diffuse_and_map(pdb = single_pdb_name,
                                supercell_num = 2,
                                size_h = p.params.size_h,
                                size_k = p.params.size_k,
                                size_l = p.params.size_l,
                                Ncpu = p.params.Ncpu,
                                write_pdb = True,
                                l = l)
        mtz_name = 'single_sc.mtz'
        pdb_name = 'single_sc.pdb'
        os.system('rm supercell_out_1.pdb')
        os.system('mv supercell_out_0.pdb {}'.format(pdb_name))
        os.system('mv single_pdb.mtz {}'.format(mtz_name))
        os.system('mv single_pdb_IDFF.map single_sc_IDFF.map')
        mtz_list.append(mtz_name)
        pdb_list.append(pdb_name)
        
###########################################################################
#            Calculate perfect scattering from single_pdb in P1           #
###########################################################################

    if p.params.single_P1:
        l.process_message('single_pdb P1 operation...')
        # use EXP class to calculate diffuse scattering
        ex = EXPClass()
        ex.calc_diffuse_and_map(pdb = single_pdb_name,
                                supercell_num = 2,
                                size_h = 1,
                                size_k = 1,
                                size_l = 1,
                                Ncpu = p.params.Ncpu,
                                write_pdb = True,
                                l = l)
        mtz_name = 'single_P1.mtz'
        pdb_name = 'single_P1.pdb'
        os.system('rm supercell_out_1.pdb')
        os.system('mv supercell_out_0.pdb {}'.format(pdb_name))
        os.system('mv single_pdb.mtz {}'.format(mtz_name))
        os.system('mv single_pdb_IDFF.map single_P1_IDFF.map')
        mtz_list.append(mtz_name)
        pdb_list.append(pdb_name)
        
###########################################################################
#                Calculate diffuse scattering from rb_pdb                 #
###########################################################################

    if p.params.rb_sc:
        l.process_message('rb_pdb supercell operation...')
        # use EXP class to calculate diffuse scattering
        ex = EXPClass()
        ex.calc_diffuse_and_map(pdb = p.output.rb_pdb_out,
                                supercell_num = 100,
                                size_h = p.params.size_h,
                                size_k = p.params.size_k,
                                size_l = p.params.size_l,
                                Ncpu = p.params.Ncpu,
                                write_pdb = True,
                                l = l)
        mtz_name = 'rb_sc.mtz'
        pdb_name = 'rb_sc.pdb'
        os.system('mv supercell_out_0.pdb {}'.format(pdb_name))
        os.system('rm supercell_out_*.pdb')
        os.system('mv rb_pdb_out.mtz {}'.format(mtz_name))
        os.system('mv rb_pdb_out_IDFF.map rb_sc_IDFF.map')
        mtz_list.append(mtz_name)
        pdb_list.append(pdb_name)
        
###########################################################################
#            Calculate diffuse scattering from rb_pdb in P1               #
###########################################################################

    if p.params.rb_P1:
        l.process_message('rb_pdb P1 operation...')
        # use EXP class to calculate diffuse scattering
        ex = EXPClass()
        ex.calc_diffuse_and_map(pdb = p.output.rb_pdb_out,
                                supercell_num = 100,
                                size_h = 1,
                                size_k = 1,
                                size_l = 1,
                                Ncpu = p.params.Ncpu,
                                write_pdb = True,
                                l = l)
        mtz_name = 'rb_P1.mtz'
        pdb_name = 'rb_P1.pdb'        
        os.system('mv supercell_out_0.pdb {}'.format(pdb_name))
        os.system('rm supercell_out_*.pdb')
        os.system('mv rb_pdb_out.mtz {}'.format(mtz_name))
        os.system('mv rb_pdb_out_IDFF.map rb_P1_IDFF.map')
        mtz_list.append(mtz_name)
        pdb_list.append(pdb_name)
        
###########################################################################
#                      Prep intensity files for er                        #
###########################################################################

    if p.params.add_sigma:
        # Modify .mtz files by adding a SIGI column! This will be sqrt(I)
        l.process_message('Adding sigmas (sqrt(I) to mtz files...')
        if len(mtz_list) == 0:
            # Dev option for when fft steps are skipped
            mtz_list = ['single_sc.mtz',
                        'single_P1.mtz',
                        'rb_sc.mtz',
                        'rb_P1.mtz']
        for mtz_file in mtz_list:
            l.show_info('Processing mtz file: {}'.format(mtz_file))
            # Read mtz, create mtz_object
            mtz_object = mtz.object(mtz_file)
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
            # Write new MTZ to file
            mtz_dataset.mtz_object().write("{}_SIG.mtz".format(mtz_file[:-4]))

###########################################################################
#                For all PDB's add B-fact noise                           #
###########################################################################

    # Add noise (0-0.1) in B-factor column, this allows for TLS fitting without
    # influencing the ensemble refinement
    if p.params.b_fact_noise:
        l.process_message('Adding noise to B-factor column')
        # Dev-option:
        if len(pdb_list) == 1:
            pdb_list = ['single_sc.pdb',
                        'single_P1.pdb',
                        'rb_sc.pdb',
                        'rb_P1.pdb']
        # Loop over all PDB's
        for pdb_file in pdb_list:
            l.show_info('Adding noise to {}'.format(pdb_file))
            # Read PDB
            pdb_f = PDBClass(fname=pdb_file, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
            # Loop over all atoms
            for atom in pdb_f.hierarchy.models()[0].atoms():
                atom.b = random.uniform(0.0,0.1)
            # Write PDB to file
            pdb_f.write_hierarchy_to_pdb(output_hierarchy=pdb_f.hierarchy,
                                         out_name='{}_B_noise.pdb'.format(pdb_file[:-4]))
            
###########################################################################
#                Prep input for ensemble refinement                       #
###########################################################################

    # P1
    single_P1_vs_single_P1 = 'phenix.ensemble_refinement single_P1_B_noise.pdb single_P1_SIG.mtz output_file_prefix=single_P1_vs_single_P1  params'
    single_P1_vs_rb_P1 = 'phenix.ensemble_refinement single_P1_B_noise.pdb rb_P1_SIG.mtz output_file_prefix=single_P1_vs_rb_P1 params'
    rb_P1_vs_single_P1 = 'phenix.ensemble_refinement rb_P1_B_noise.pdb single_P1_SIG.mtz output_file_prefix=rb_P1_vs_rb_P1 params'
    rb_P1_vs_rb_P1 = 'phenix.ensemble_refinement rb_P1_B_noise.pdb rb_P1_SIG.mtz output_file_prefix=rb_P1_vs_single_P1 params'
    
    # SC
    single_sc_vs_single_sc = 'phenix.ensemble_refinement single_sc_B_noise.pdb single_sc_SIG.mtz output_file_prefix=single_sc_vs_single_sc params'
    single_sc_vs_rb_sc = 'phenix.ensemble_refinement single_sc_B_noise.pdb rb_sc_SIG.mtz output_file_prefix=single_sc_vs_rb_sc params'
    rb_sc_vs_single_sc = 'phenix.ensemble_refinement rb_sc_B_noise.pdb single_sc_SIG.mtz output_file_prefix=rb_sc_vs_single_sc params'
    rb_sc_vs_rb_sc = 'phenix.ensemble_refinement rb_sc_B_noise.pdb rb_sc_SIG.mtz output_file_prefix=rb_sc_vs_rb_sc params'

    # Running parameters (right now for TESTING!!!!!) adjust tx for real runs!!!
    params_commands = '''
ensemble_refinement {
  max_ptls_cycles=1
  tls_group_selections = all
  ptls = 0.0
  tx = 1.0
  equilibrium_n_tx = 2
  acquisition_block_n_tx = 4
  number_of_aquisition_periods = 5
  cartesian_dynamics.stop_cm_motion = False
  ordered_solvent_update = False
  ensemble_reduction = False
  output_running_kinetic_energy_in_occupancy_column = True
}
input.xray_data.labels = ITOT,SIGITOT
input.xray_data.r_free_flags.generate=True
 '''

    # Write parameter file
    with open('params','w') as f:
        print >> f, params_commands

###########################################################################
#                   Start simulations (parallel and screened)             #
###########################################################################

    # P1
    com = 'screen -dmSL {} {}'.format('single_P1_vs_single_P1',single_P1_vs_single_P1)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('single_P1_vs_rb_P1',single_P1_vs_rb_P1)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('rb_P1_vs_single_P1',rb_P1_vs_single_P1)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('rb_P1_vs_rb_P1',rb_P1_vs_rb_P1)
    os.system(com)
    # SC
    com = 'screen -dmSL {} {}'.format('single_sc_vs_single_sc',single_sc_vs_single_sc)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('single_sc_vs_rb_sc',single_sc_vs_rb_sc)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('rb_sc_vs_single_sc',rb_sc_vs_single_sc)
    os.system(com)
    com = 'screen -dmSL {} {}'.format('rb_sc_vs_rb_sc',rb_sc_vs_rb_sc)
    os.system(com)
    
    return l
