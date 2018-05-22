#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

from scitbx.array_family import flex

import phil as additivity_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
from scud.mtz.MTZ import mtz_correlation

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
    l.title("exp.additivity module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = additivity_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.additivity

###########################################################################
#                           Read mtz_files                                #
###########################################################################

    l.process_message('Reading mtz_1...')
    mtz_1 = MTZClass(fname = p.input.mtz_1,
                     array_name = 'IDFF')

    l.process_message('Reading mtz_2...')
    mtz_2 = MTZClass(fname = p.input.mtz_2,
                     array_name = 'IDFF')

#    l.process_message('Calculating R-value and CC...')
#    mtz_correlation(mtz_1 = mtz_1,
#                    mtz_2= mtz_2,
#                    super_1 = 1,
#                    super_2 = 1,
#                    l = l)

    temp_array = mtz_1.miller_array.deep_copy()
    t = temp_array.data() - mtz_2.miller_array.data()
    temp_miller = mtz_1.miller_array.customized_copy(data=t)
 
    mtz_1.write_array_as_mtz(array = temp_miller,
                             label = 'IDFF',
                             mtz_out_name=p.output.mtz_out)
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    return l
    
def old():
    print_title('tst_additivity')
    single = sys.argv[1]
    ensemble = sys.argv[2]
    #### INIT PARAMS
    rotation_sigma_deg = 3.5
    mu = 0
    tst_B = 15
    translation_sigma = np.sqrt(tst_B / (8 * (np.pi**2) ))
    ### Read PDB's
    pdb_single = PDBClass(fname=single, selection="not (resname HOH) and not (resname CL) and not (resname NA) and not (altid B) and not (altid C)")
#    pdb_ensemble = PDBClass(fname=ensemble, selection="(altid " " or altid A) and not (resname HOH) and not (resname CL) and not (resname NA)")
    pdb_ensemble = PDBClass(fname=ensemble, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    #### Set ensemble size to be equal to # models in input ensemble:
    ensemble_size = int(pdb_ensemble.hierarchy.models()[-1].id)
    #### RMSF in B-Factors of original ensemble
    pdb_ensemble.start_or_attach_hierarchy()
    for i in range(ensemble_size):
        mdl = pdb_ensemble.detached_model_copy(n=i)
        pdb_ensemble.start_or_attach_hierarchy(detached_model=mdl)
    out_name_ens_plot = 'ensembles/'+datetime_string()+'_ensemble_simple.json'
    calc_RMSF_and_write(pdb_ensemble,lab="Ensemble Simple",fname=out_name_ens_plot,col=(0.2,0.9,0.2))
    #### COM for rotations:
    center_of_mass =  pdb_single.calculate_center_of_mass(pdb_single.hierarchy)
    #### Generate trans_list and rot_list
    translation_vector_list = pdb_single.create_translation_vector_list(mu=mu,
                                                                    translation_sigma=translation_sigma,
                                                                        ensemble_size=ensemble_size)
    rotation_list = pdb_single.create_rotation_matrix_list(rotation_sigma_deg=rotation_sigma_deg,
                                                           ensemble_size=ensemble_size)
    #### Set all B-Values to 0, B is put in motion of atoms
    pdb_single.set_B_zero()
    
    print '#### Single Translation ####'
    #### Initialize, by generating a new/empty heirarchy
    pdb_single.start_or_attach_hierarchy()
    for translation_vector in translation_vector_list:
        detached_model = pdb_single.detached_model_copy()        
        translated_model = pdb_single.translate_model(detached_model,
                                                      translation_vector)
        pdb_single.start_or_attach_hierarchy(detached_model=translated_model)
    #### Write PDB to file
    out_name_trans_pdb = 'ensembles/'+datetime_string()+'_single_translation.pdb'
    trans_plot_name = out_name_trans_pdb[:-4]+".json"
    calc_RMSF_and_write(pdb_single,lab="Single Translation",fname=trans_plot_name,col=(0.2,0.4,0.4))
    pdb_single.write_hierarchy_to_pdb(out_name=out_name_trans_pdb)

    print '#### Single Rotation ####'
    pdb_single.start_or_attach_hierarchy()
    for rotation_matrix in rotation_list:
        detached_model = pdb_single.detached_model_copy()        
        rotated_model = pdb_single.rotate_model(detached_model,
                                                rotation_matrix,
                                                center_of_mass)
        pdb_single.start_or_attach_hierarchy(detached_model=rotated_model)
    #### Write PDB to file
    out_name_rot_pdb = 'ensembles/'+datetime_string()+'_single_rotation.pdb'
    rot_plot_name = out_name_rot_pdb[:-4]+".json"
    calc_RMSF_and_write(pdb_single,lab="Single Rotation",fname=rot_plot_name,col=(0.9,0.9,0.2))
    pdb_single.write_hierarchy_to_pdb(out_name=out_name_rot_pdb)
    
    print '#### Single Mix ####'
    pdb_single.start_or_attach_hierarchy()
    for rotation_matrix,translation_vector in zip(rotation_list,translation_vector_list):
        detached_model = pdb_single.detached_model_copy()
        rotated_model = pdb_single.rotate_model(detached_model,
                                                rotation_matrix,
                                                center_of_mass)
        translated_model = pdb_single.translate_model(rotated_model,
                                                      translation_vector)
        pdb_single.start_or_attach_hierarchy(detached_model=translated_model)
    #### Write PDB to file
    out_name_mix_pdb = 'ensembles/'+datetime_string()+'_single_mix.pdb'
    mix_plot_name = out_name_mix_pdb[:-4]+".json"
    calc_RMSF_and_write(pdb_single,lab="Single Mix",fname=mix_plot_name,col=(0.2,0.9,0.9))
    pdb_single.write_hierarchy_to_pdb(out_name=out_name_mix_pdb)

    print '#### Mix ensemble #####'
    # Quite different, need to extract a hundred models from the ensemble (randomly choose)
    # And apply rot/trans operations on them
    pdb_ensemble.set_B_zero()
    
    pdb_ensemble.start_or_attach_hierarchy()
    n = 0
    for rotation_matrix,translation_vector in zip(rotation_list,translation_vector_list):
        detached_model = pdb_ensemble.detached_model_copy(n=n)
        rotated_model = pdb_ensemble.rotate_model(detached_model,
                                                rotation_matrix,
                                                center_of_mass)
        translated_model = pdb_ensemble.translate_model(rotated_model,
                                                        translation_vector)
        pdb_ensemble.start_or_attach_hierarchy(detached_model=translated_model)
        n+=1
    #### Write PDB to file
    out_name_ens_mix_pdb = 'ensembles/'+datetime_string()+'_ensemble_mix.pdb'
    ens_mix_plot_name = out_name_ens_mix_pdb[:-4]+".json"
    calc_RMSF_and_write(pdb_ensemble,lab="Ensemble Mix",fname=ens_mix_plot_name,col=(0.2,0.2,0.9))
    pdb_ensemble.write_hierarchy_to_pdb(out_name=out_name_ens_mix_pdb)

    #### Plot all RMSF as B-Factors
    plot_all_lines([out_name_ens_plot,trans_plot_name,rot_plot_name,mix_plot_name,ens_mix_plot_name])
    
def plot_all_lines(fnames):
    import json
    dat_list = []
    for f in fnames:
        with open(f) as fl:
            dat = json.load(fl)
            dat_list.append(dat)
    p = Plot()
    p.line_data = dat_list
    p.lines_plot()
            
def calc_RMSF_and_write(pdb_obj,lab=None,fname=None,col=(0,0,0)):
    '''
    Calculate ensemble RMSF as B-Factors and write to file
    '''
    x_dat, y_dat = pdb_obj.ensemble_to_B_factor(ensemble_hierarchy=pdb_obj.new_hierarchy)
    p = Plot()
    p.collect_line_data(xlabel='Residue Number',
                        ylabel='B-factor ($\AA^{2}$)',
                        label=lab,
                        color=col,
                        x=x_dat,
                        y=y_dat,
                        fname=fname)
    p.write_line_data()


    
# Run Program
#run()
