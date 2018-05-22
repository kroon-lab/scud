#!/usr/bin/env cctbx.python

import numpy as np
import sys,os

from scitbx.array_family import flex

from scud.general.method_library import print_title
from scud.general.method_library import datetime_string
from scud.pdb.PDB import PDBClass
from scud.pdb.PDB import Optimiser
from scud.pdb.PDB import CA_Optimiser
from scud.general.Plot import Plot
import phil as rigid_body_phil
from libtbx.phil import parse
from scud.general.log_writer import Log

def run(args=None, l=None):
    """
    Generate an ensemble describing types of rigid body motion,
    input can be single structure or ensemble describng 'internal' motion

    # TODO:
    * Add alternative rigid_body method where each motion will have the same
      average B-factor as the input model (measured model)

    * Add cutoff for too high B-factors in calculation average B (ensemble 
      refinement can grosely overestimate B-factor in RMSF)
    """
###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    if l == None:
        l = Log(log_name='rigid_body_log.txt')
    l.title("pdb.rigid_body module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    working_params = rigid_body_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.rigid_body

    l.process_message('rb_types: er, trans, rot, mix')

    # er type, fit sigma rot and trans to btls in B-factor column of ensemble
    # Apply these to the input ensemble to obtain an ensemble describing both
    # Atomic motion and TLS motions
    if p.params.rb_type == 'er':
        l.process_message('rb_type = er')
        # Use for fitting trans and rot sigma:
        pdb_single = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Extract btls from and apply final rot and trans on
        btls_pdb = PDBClass(fname=p.input.ensemble_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Plotting
        classic_B = [at.b for at in pdb_single.hierarchy.atoms() if at.name == ' CA ']
        # Plottign x-axis
        res_num = [int(resid.resseq) for resid in pdb_single.hierarchy.models()[0].chains()[0].residues()]
        # Plotting class
        pt = Plot()
        # Plot classic B-factors as reference
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='classic B-factors',
                             color=(0.2,0.2,0.8),
                             x=res_num,
                             y=classic_B,
                             fname="er_"+p.output.plot_out[:-4]+"_classic_B.json")
        # BTLS, used for plotting and target function for fitten
        btls = [at.b for at in btls_pdb.hierarchy.models()[0].atoms() if at.name == ' CA ']
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='BTLS',
                             color=(0.8,0.2,0.8),
                             x=res_num,
                             y=btls,
                             fname="er_"+p.output.plot_out[:-4]+"_btls_B.json")
        # Write json with line data
        pt.write_line_data()
        # Plot line data
        pt.lines_plot(plotfile_name='er_start.eps')
        l.process_message('Perform actual operations...')
        # Assign target b-factor trace (CA)
        pdb_single.target_B_CA = np.array(btls)
        # Set all B-Values to 0, B is put in motion of atoms
        pdb_single.set_B_zero()
        # Set occupancy to 1
        pdb_single.set_occ()
        # COM for rotations
        center_of_mass =  pdb_single.calculate_center_of_mass()
        # For creating multi models
        pdb_single.model_list = pdb_single.create_random_model_list(l = l,
                                                                    ensemble_size = p.input.ensemble_size)
        # Call simplex class
        start_trans_sigma, start_rot_sigma, = 0.5, 0.4
        simp_opt = CA_Optimiser(trans_sigma = start_trans_sigma,
                                rot_sigma = start_rot_sigma,
                                pdb_object = pdb_single,
                                center_of_mass = center_of_mass,
                                rb_type = 'mix',
                                l = l)    
        ls = list(simp_opt.result)
        # Final results
        trans_sigma, rot_sigma = p.params.scale_t*ls[0],p.params.scale_r*ls[1]
        # Apply to input esnemble
        btls_pdb.model_list = btls_pdb.create_random_model_list(l = l,
                                                                ensemble_size = p.input.ensemble_size)
        btls_pdb.set_B_zero()
        # Set occupancy to 1
        btls_pdb.set_occ()
        center_of_mass =  btls_pdb.calculate_center_of_mass()
        # Create final ensemble, ER result + added btls motion
        btls_pdb.make_rot_trans_ens(trans_only = False,
                                    rot_only = False,
                                    trans_sigma = trans_sigma,
                                    rot_sigma = rot_sigma,
                                    center_of_mass = center_of_mass,
                                    l = l)
        # Write final ensemble
        btls_pdb.write_hierarchy_to_pdb(out_name='er_ens_rb.pdb')
        # Calculate B-factor from ensemble rmsf
        atom,b_total = btls_pdb.ensemble_to_B_factor(ensemble_hierarchy=btls_pdb.new_hierarchy,
                                                     ca_only = True)
        # Plot B-total/B-new
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='rb-ens',
                             color=(0.8,0.8,0.2),
                             x=res_num,
                             y=b_total,
                             fname="er_"+p.output.plot_out[:-4]+"_new_B.json")
        # Plot final results
        pt.write_line_data()
        pt.lines_plot(plotfile_name='er_finish.eps')

    # trans, fit translation only to B-factor trace of input pdb
    if p.params.rb_type == 'trans':
        l.process_message('rb_type = trans')
        # Use for fitting trans and rot sigma:
        pdb_single = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Plot 
        classic_B = [at.b for at in pdb_single.hierarchy.atoms() if at.name == ' CA ']
        # Plottign x-axis
        res_num = [int(resid.resseq) for resid in pdb_single.hierarchy.models()[0].chains()[0].residues()]
        # Plotting class
        pt = Plot()
        # Plot classic B-factors as reference
        pdb_single.target_B_CA = np.array(classic_B)
        # Set all B-Values to 0, B is put in motion of atoms
        pdb_single.set_B_zero()
        # Set occupancy to 1
        pdb_single.set_occ()
        # COM for rotations
        center_of_mass =  pdb_single.calculate_center_of_mass()
        # For creating multi models
        pdb_single.model_list = pdb_single.create_random_model_list(l = l,
                                                                    ensemble_size = p.input.ensemble_size)
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='classic B-factors',
                             color=(0.2,0.2,0.8),
                             x=res_num,
                             y=classic_B,
                             fname="trans_"+p.output.plot_out[:-4]+"_classic_B.json")
        start_trans_sigma, start_rot_sigma, = 0.5, 0.0
        simp_opt = CA_Optimiser(trans_sigma = start_trans_sigma,
                                rot_sigma = start_rot_sigma,
                                pdb_object = pdb_single,
                                center_of_mass = center_of_mass,
                                rb_type = 'trans',
                                l = l)
        ls = list(simp_opt.result)
        # Final results
        trans_sigma, rot_sigma = p.params.scale_t*ls[0],p.params.scale_r*ls[1]
        pdb_single.make_rot_trans_ens(trans_only = True,
                                      rot_only = False,
                                      trans_sigma = trans_sigma,
                                      rot_sigma = rot_sigma,
                                      center_of_mass = center_of_mass,
                                      l = l)
        pdb_single.write_hierarchy_to_pdb(out_name='trans_rb.pdb')
        # Calculate B-factor from ensemble rmsf
        atom,b_total = pdb_single.ensemble_to_B_factor(ensemble_hierarchy=pdb_single.new_hierarchy,
                                                       ca_only = True)
        # Plot B-total/B-new
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='rb-trans',
                             color=(0.2,0.2,0.8),
                             x=res_num,
                             y=b_total,
                             fname="trans_"+p.output.plot_out[:-4]+"_new_B.json")
        # Plot final results
        pt.write_line_data()
        pt.lines_plot(plotfile_name='trans_finish.eps')
        # trans, fit translation only to B-factor trace of input pdb
    if p.params.rb_type == 'rot':
        l.process_message('rb_type = rot')
        # Use for fitting trans and rot sigma:
        pdb_single = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Plot 
        classic_B = [at.b for at in pdb_single.hierarchy.atoms() if at.name == ' CA ']
        # Plottign x-axis
        res_num = [int(resid.resseq) for resid in pdb_single.hierarchy.models()[0].chains()[0].residues()]
        # Plotting class
        pt = Plot()
        # Plot classic B-factors as reference
        pdb_single.target_B_CA = np.array(classic_B)
        # Set all B-Values to 0, B is put in motion of atoms
        pdb_single.set_B_zero()
        # Set occupancy to 1
        pdb_single.set_occ()
        # COM for rotations
        center_of_mass =  pdb_single.calculate_center_of_mass()
        # For creating multi models
        pdb_single.model_list = pdb_single.create_random_model_list(l = l,
                                                                    ensemble_size = p.input.ensemble_size)
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='classic B-factors',
                             color=(0.2,0.8,0.2),
                             x=res_num,
                             y=classic_B,
                             fname="rot_"+p.output.plot_out[:-4]+"_classic_B.json")
        start_trans_sigma, start_rot_sigma, = 0.0, 0.4
        simp_opt = CA_Optimiser(trans_sigma = start_trans_sigma,
                                rot_sigma = start_rot_sigma,
                                pdb_object = pdb_single,
                                center_of_mass = center_of_mass,
                                rb_type = 'rot',
                                l = l)
        ls = list(simp_opt.result)
        # Final results
        trans_sigma, rot_sigma = p.params.scale_t*ls[0],p.params.scale_r*ls[1]
        pdb_single.make_rot_trans_ens(trans_only = False,
                                      rot_only = True,
                                      trans_sigma = trans_sigma,
                                      rot_sigma = rot_sigma,
                                      center_of_mass = center_of_mass,
                                      l = l)
        pdb_single.write_hierarchy_to_pdb(out_name='rot_rb.pdb')
        # Calculate B-factor from ensemble rmsf
        atom,b_total = pdb_single.ensemble_to_B_factor(ensemble_hierarchy=pdb_single.new_hierarchy,
                                                       ca_only = True)
        # Plot B-total/B-new
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='rb-rot',
                             color=(0.8,0.2,0.2),
                             x=res_num,
                             y=b_total,
                             fname="rot_"+p.output.plot_out[:-4]+"_new_B.json")
        # Plot final results
        pt.write_line_data()
        pt.lines_plot(plotfile_name='rot_finish.eps')
    if p.params.rb_type == 'mix':
        l.process_message('rb_type = mix')
        # Use for fitting trans and rot sigma:
        pdb_single = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        # Plot 
        classic_B = [at.b for at in pdb_single.hierarchy.atoms() if at.name == ' CA ']
        # Plottign x-axis
        res_num = [int(resid.resseq) for resid in pdb_single.hierarchy.models()[0].chains()[0].residues()]
        # Plotting class
        pt = Plot()
        # Plot classic B-factors as reference
        pdb_single.target_B_CA = np.array(classic_B)
        # Set all B-Values to 0, B is put in motion of atoms
        pdb_single.set_B_zero()
        # Set occupancy to 1
        pdb_single.set_occ()
        # COM for rotations
        center_of_mass =  pdb_single.calculate_center_of_mass()
        # For creating multi models
        pdb_single.model_list = pdb_single.create_random_model_list(l = l,
                                                                    ensemble_size = p.input.ensemble_size)
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='classic B-factors',
                             color=(0.2,0.2,0.8),
                             x=res_num,
                             y=classic_B,
                             fname="mix_"+p.output.plot_out[:-4]+"_classic_B.json")
        start_trans_sigma, start_rot_sigma, = 0.1, 0.4
        simp_opt = CA_Optimiser(trans_sigma = start_trans_sigma,
                                rot_sigma = start_rot_sigma,
                                pdb_object = pdb_single,
                                center_of_mass = center_of_mass,
                                rb_type = 'mix',
                                l = l)
        ls = list(simp_opt.result)
        # Final results
        trans_sigma, rot_sigma = p.params.scale_t*ls[0],p.params.scale_r*ls[1]
        pdb_single.make_rot_trans_ens(trans_only = False,
                                      rot_only = False,
                                      trans_sigma = trans_sigma,
                                      rot_sigma = rot_sigma,
                                      center_of_mass = center_of_mass,
                                      l = l)
        pdb_single.write_hierarchy_to_pdb(out_name='mix_rb.pdb')
        # Calculate B-factor from ensemble rmsf
        atom,b_total = pdb_single.ensemble_to_B_factor(ensemble_hierarchy=pdb_single.new_hierarchy,
                                                       ca_only = True)
        # Plot B-total/B-new
        pt.collect_line_data(xlabel='CA',
                             ylabel='B-factor ($\AA^{2}$)',
                             label='rb-mix',
                             color=(0.2,0.2,0.2),
                             x=res_num,
                             y=b_total,
                             fname="mix_"+p.output.plot_out[:-4]+"_new_B.json")
        # Plot final results
        pt.write_line_data()
        pt.lines_plot(plotfile_name='mix_finish.eps')
    return l
        
def old_not_old():
###########################################################################
#                      Read and analyse pdb_in                            #
###########################################################################

    # p.input.b_ref !!!!!
    l.process_message('Reading and analyzing PDB file / Structure...')
    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    # Model list is list with model numbers randomly sorted
    model_list = pdb_object.create_random_model_list(l = l,
                                                     ensemble_size = p.input.ensemble_size)
    l.show_info('Number of models in PDB: {}'.format(pdb_object.model_num))
    # Reference B-factor for fitting rigid body motion to
    if p.input.b_ref != None:
        # Checking B-Factors, will be used to determine translational magnitude
        l.process_message('Using reference PDB to determine B-facotr profile to fit to...')
        pdb_ref = PDBClass(fname=p.input.b_ref, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        target_b_factor = pdb_ref.average_B(make_mask = True)
    if p.input.b_ext != None:
        l.process_message('Using external B-factor: {}...'.format(p.input.b_ext))
        target_b_factor = p.input.b_ext
    else:
        target_b_factor = pdb_object.average_B(make_mask = True)
    l.show_info('Target B-Factor (A^2): {:.2f}'.format(target_b_factor))
    x_ca_b_factor,y_ca_b_factor = pdb_object.extract_ca_b_factors(pdb_object.hierarchy.deep_copy())
    # Calculate RMSF which will be used to generate the translation list
    translation_sigma = pdb_object.B_factor_to_RMSF(b = target_b_factor) 
    l.show_info('translation_sigma (A): {:.2f}'.format(translation_sigma))
    l.show_info('rotation_sigma (degrees): {:.2f}'.format(p.params.start_rotation_sigma))
    center_of_mass =  pdb_object.calculate_center_of_mass()
    # Center of mass needed for rotational operations
    l.show_info('center of mass (A): {:.2f} {:.2f} {:.2f}'.format(*center_of_mass))
    

###########################################################################
#                      Start translation and rotation                     #
###########################################################################

    l.process_message('Perform actual operations...')
    # Set all B-Values to 0, B is put in motion of atoms
    pdb_object.set_B_zero()
    # Set occupancy to 1
    pdb_object.set_occ()
    pdb_object.model_list = model_list
    if p.input.b_ref != None:
        pdb_object.mask = pdb_ref.mask
        pdb_object.masked_target_b = pdb_ref.masked_target_b
    if p.input.b_ext != None:
        pdb_object.mask = np.array([True for i in pdb_object.hierarchy.models()[0].atoms()])
        pdb_object.masked_target_b = np.array([target_b_factor
                                               for i in pdb_object.hierarchy.models()[0].atoms()])
        
###########################################################################
#                Simplex to minimum rms B-factor profile                  #
###########################################################################

    trans_sigma, rot_sigma = translation_sigma, p.params.start_rotation_sigma
    simp_opt = Optimiser(trans_sigma = trans_sigma,
                         rot_sigma = rot_sigma,
                         trans_only = p.params.translation_only,
                         rot_only = p.params.rotation_only,
                         pdb_object = pdb_object,
                         center_of_mass = center_of_mass,
                         l = l)    
    # Initialize, by generating a new/empty heirarchy
    ls = list(simp_opt.result)
    trans_sigma, rot_sigma = ls[0],ls[1]
    #### Final rotation/translation
    pdb_object.make_rot_trans_ens(trans_only = p.params.translation_only,
                                  rot_only = p.params.rotation_only,
                                  trans_sigma = trans_sigma,
                                  rot_sigma = rot_sigma,
                                  center_of_mass = center_of_mass,
                                  l = l)
    atom, b = pdb_object.ensemble_to_B_factor(ensemble_hierarchy=pdb_object.new_hierarchy,
                                             ca_only = False)
    b_in = [a.b for a in pdb_object.hierarchy.atoms()]
    pt = Plot()
    pt.collect_line_data(xlabel='Atom Number',
                         ylabel='B-factor ($\AA^{2}$)',
                         label='rmsf B-factors',
                         color=(0.2,0.2,0.8),
                         x=atom,
                         y=b,
                         fname="tmp_b_from_ens.json")
    pt.collect_line_data(xlabel='Atom Number',
                         ylabel='B-factor ($\AA^{2}$)',
                         label='Input B-factors',
                         color=(0.2,0.2,0.2),
                         x=range(len(b_in)),
                         y=b_in,
                         fname="tmp_b_in.json")
    pt.write_line_data()
    l.show_info('plot file name: {}'.format('tst.eps'))
    pt.lines_plot(plotfile_name='tst.eps')
    
###########################################################################
#           Finishing up, writing output files/plots                      #
###########################################################################

    l.process_message('Writing new ensemble to file...')
    l.show_info('file name: {}'.format(p.output.pdb_out))
    #### Write PDB to file
    pdb_object.write_hierarchy_to_pdb(out_name=p.output.pdb_out)

###########################################################################
#                       Analysing new ensemble                            #
###########################################################################

    l.process_message('Analyzing new ensemble...')
    #### Calculate B factor from ensemble RMSF
    x_ens_b,y_ens_b = pdb_object.ensemble_to_B_factor(ensemble_hierarchy=pdb_object.start_or_attach_hierarchy(return_hierarchy=True))
    l.process_message('Plotting B-Factors...')
    #### Plot
    pt = Plot()
    pt.collect_line_data(xlabel='Residue Number',
                         ylabel='B-factor ($\AA^{2}$)',
                         label='classic B-factors',
                         color=(0.2,0.2,0.8),
                         x=x_ca_b_factor,
                         y=y_ca_b_factor,
                         fname=p.output.plot_out[:-4]+"_classic_B.json")
    pt.collect_line_data(xlabel='Residue Number',
                         ylabel='B-factor ($\AA^{2}$)',
                         label='Ensemble B (RMSF)',
                         color=(0.8,0.4,0.0),
                         x=x_ens_b,
                         y=y_ens_b,
                         fname=p.output.plot_out[:-4]+"_ensemble_B.json")
    pt.write_line_data()
    l.show_info('plot file name: {}'.format(p.output.plot_out))
    pt.lines_plot(plotfile_name=p.output.plot_out)
    
    return l
