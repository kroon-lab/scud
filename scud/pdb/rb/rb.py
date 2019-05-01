#!/usr/bin/env cctbx.python

import numpy as np
import sys,os
import matplotlib.pyplot as plt
import matplotlib

from scitbx.array_family import flex

import phil as rb_phil
from libtbx.phil import parse
from scud.general.log_writer import Log

from rb_pdb import RB_PDB
from rb_pdb import RB_Optimiser
from rb_pdb import RB_Aniso_Optimiser

def run(args = None,
        l = None):

    '''
    - Extract B-factor profile from target_pdb
    - Create template for enesmeble generation
    - Simplex minimize translation and rotation parameters with target B
      as a target

    Classes:
       * RB_PDB, PDB reading and manipulations

    Input:
       * target_pdb, B-factors will be used to fit rigid body motion against
       * template_pdb, will be used to generate rigid body ensemble, may be
                       single or multi model pdb file
       * rb_type, trans/mix/rot

    Parameters:
       * filter_target_b, filter high or low B-factors before fitted
       * filter_n_sighma, number of std's before B-factor is filtered
       * ensemble_size, size of final ensemble

    Output:
      * Plot with target and fitted B-factor traces
      * Plot showing distribution of trans/rot vectors/angles
      * PDB file containing new ensemble


    TODO:
      * Add 'new' translation algorithm

    '''

###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    if l == None:
        l = Log(log_name='rb_log.txt')
    l.title("pdb.rb module")

###########################################################################
#                      Process input Params                               #
###########################################################################

    # Read input parameters
    working_params = rb_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.rb
    
###########################################################################
#                             Set Seed                                    #
###########################################################################

    # Set seed for random number generators used for translation and rotation
    np.random.seed(p.params.seed)


###########################################################################
#                          Read PDB file                                  #
###########################################################################

    ####  Read PDB ####

    # Check if a different target pdb is supplies
    if p.input.target_pdb == None:
        target_pdb_fname = p.input.template_pdb
        
    else:
        target_pdb_fname = p.input.target_pdb

    l.process_message('Reading pdb file used for B-factor fitting')

    # (extract (selected) atoms and symmetry, create hierachy)
    target_pdb = RB_PDB(fname = target_pdb_fname)
    target_pdb.read_pdb(selection="not (resname HOH) and not (resname CL) and not (resname NA)",
                        l = l)

    #### read CA B-factors ####

    target_pdb.extract_target_b()
    l.show_info('Average CA B-factor: {}'
                .format(np.mean(target_pdb.target_b)))

    #### Obtain residue numbers ####

    # For plotting later on:
    residue_numbers = [int(resid.resseq) for resid in target_pdb.hierarchy.models()[0].chains()[0].residues()]
    
    

###########################################################################
#                          Prep B-factors                                 #
###########################################################################

    #### Mask b bactors ###

    l.process_message('Prepping B-factors for simplex fitting procedure...')

    #(based on deviation from median)

    if p.params.filter_target_b:

        # Get sigma from input
        n_sigma = p.params.filter_n_sigma

        # Print info
        l.show_info('n_sigma for B-factor filtering: {}'.format(n_sigma))

        # Create B-factor mask
        target_pdb.create_b_factor_mask(n_sigma=n_sigma)
        filter_b = target_pdb.target_b[target_pdb.mask]
#        print 'mask',target_pdb.mask
        print 'filter_b', filter_b

        # Print info
        l.show_info('Average CA B-factor after masking: {}'
                    .format(np.mean(filter_b)))
        l.show_info('B-factor of {} converts to RMSF of {}'.format(np.mean(filter_b),
                                                                   b_to_rmsf(np.mean(filter_b))))
    # If no filtering is needed create 'True mask'
    if not p.params.filter_target_b:

        # Create list of ones of len(target_b) convert to type boolean
        target_pdb.mask = np.ones(len(target_pdb.target_b)).astype(bool)

###########################################################################
#                          Prep template_PDB                              #
###########################################################################

    #### Read template PDB ###

    l.process_message('Preppign template PDB...')

    template_pdb = RB_PDB(fname = p.input.template_pdb)
    template_pdb.read_pdb(selection="not (resname HOH) and not (resname CL) and not (resname NA)",
                          l = l)

    #### Check for multi model ####

    # Number of models
    template_pdb.n_models = len(template_pdb.hierarchy.models())

    # If single model or less models then in input:
    if template_pdb.n_models == 1 or template_pdb.n_models < p.params.ensemble_size:
        ens_size = p.params.ensemble_size
        l.show_info('Final ensemble wil contain {} models'.format(ens_size))

    # If template contains more models then p.params.ensemble_size
    if template_pdb.n_models > p.params.ensemble_size:
        ens_size = template_pdb.n_models
        l.show_info('Final ensemble wil contain {} models'.format(ens_size))

    #### Prep template ###

    # Set B=0.0, occ=1.0
    template_pdb.set_b_factor()
    template_pdb.set_occupancy()

    #### Calculate COM ####

    template_pdb.calculate_center_of_mass()
    l.show_info('Center of mass of template PDB: {}'.format(template_pdb.com))

###########################################################################
#                           Start simplex                                 #
###########################################################################

    #### If not aniso optimize 1 trans and 1 rot parameter: ####

    if p.params.aniso == False:
        if p.input.rb_type == 'rot':

            #### Initialize start values for rb_type rot ####

            start_trans_sigma = 0.0
            l.show_info('Start trans_sigma = {}'.format(start_trans_sigma))
            start_rot_sigma = 3.

        elif p.input.rb_type == 'trans':

            #### Initialize start values for rb_type trans ####

            #### rmsf^2 = 2*sigma^2  in normal distribution ####
            
            start_trans_sigma = b_to_rmsf(np.mean(filter_b))*1.4
            l.show_info('Start trans_sigma = {}'.format(start_trans_sigma))
            start_rot_sigma = 3.

        elif p.input.rb_type == 'mix':

            #### Initialize start values for rb_type mix ####

            # Mix works better if the rotation angle is fitted once before both
            # parameters are optimize simultaniously

            l.process_message('Starting 2 step rb fitting, finding initial rotation angle...')
            start_trans_sigma = 0.4
            l.show_info('Start trans_sigma = {}'.format(start_trans_sigma))
            start_rot_sigma = 3.

            #### Optimize rotation only ####

            rb_optimizer = RB_Optimiser(trans_sigma = start_trans_sigma,
                                        rot_sigma = start_rot_sigma,
                                        template_pdb = template_pdb,
                                        ens_size = ens_size,
                                        target_b = target_pdb.target_b[target_pdb.mask],
                                        mask = target_pdb.mask,
                                        rb_type = 'rot',
                                        l = l)

            # For dual optimization use fractions of average B and optimized rotation only

            start_trans_sigma = b_to_rmsf(np.mean(filter_b))
            start_rot_sigma = 0.8 * list(rb_optimizer.result)[0]

        else:
            raise Exception('Not an rb_type')

        #### Perform simplex ####

        rb_optimizer = RB_Optimiser(trans_sigma = start_trans_sigma,
                                    rot_sigma = start_rot_sigma,
                                    template_pdb = template_pdb,
                                    ens_size = ens_size,
                                    target_b = target_pdb.target_b[target_pdb.mask],
                                    mask = target_pdb.mask,
                                    rb_type = p.input.rb_type,
                                    l = l)

        l.process_message('Simplex minimization finished...')

        #### Finish up simplex ####

        # Show simplex results
        r = list(rb_optimizer.result)
        if len(r) == 2:
            l.show_info('Results: \n\ttranslation (A): {:.2f} \n\tangle (deg): {:.2f}'.format(r[0], r[1]))
        elif len(r) == 1 and p.input.rb_type == 'rot':
            l.show_info('Results: \n\tangle (deg): {:.2f}'.format(r[0]))
        elif len(r) == 1 and p.input.rb_type == 'trans':
            l.show_info('Results: \n\ttranslation (A): {:.2f}'.format(r[0]))
        else:
            raise Exception('Invalid Value')

    #### Anisotropic fitting of rigid body motion ####

    elif p.params.aniso == True:

        #### Rational starting values of simplex ####

        # Translation start
        start_trans = b_to_rmsf(np.min(filter_b))
        # rotation start
        rot_sigma_array = np.array([2.,2.,2.])
        # center of mass start
        com = template_pdb.com

        #### Simplex Control ####

        # Step size change within simplex method
        step_change_list = [1., .75, 0.5, 0.1]

        # Start values of simple (user input)
        start_step_list = [p.params.trans_step,
                           p.params.rot_x_step,
                           p.params.rot_y_step,
                           p.params.rot_z_step,
                           p.params.com_x_step,
                           p.params.com_y_step,
                           p.params.com_z_step]
        # Translation can not be fitted anisotropically
        if p.input.rb_type == 'trans':
            l.process_message('Translation cannot be fitted anisotropically...\n quiting.')
            quit()
        
        #### Simplex rot ####

        elif p.input.rb_type == 'rot':

            #### Collect start values ####

            rot = rot_sigma_array
            com = template_pdb.com

            #### Loop over step size list, update start values ####

            for step in step_change_list:

                # Start simplex
                rb_optimizer = RB_Aniso_Optimiser(trans_sigma = 0.0,
                                                  rot_sigma_array = rot,
                                                  center_of_mass = com,
                                                  template_pdb = template_pdb,
                                                  ens_size = ens_size,
                                                  target_b = target_pdb.target_b[target_pdb.mask],
                                                  mask = target_pdb.mask,
                                                  rb_type = 'rot',
                                                  step = step,
                                                  start_step_list = start_step_list,
                                                  l = l)
                r = list(rb_optimizer.result)
                rot = r[0:3]
                com = r[3:6]

        #### Simplex mix ####

        elif p.input.rb_type == 'mix':
            
            #### Collect start values ####

            trans = start_trans
            rot = rot_sigma_array
            com = template_pdb.com

            #### Loop over step size list, update start values ####

            for step in step_change_list:

                rb_optimizer = RB_Aniso_Optimiser(trans_sigma = trans,
                                                  rot_sigma_array = rot,
                                                  center_of_mass = com,
                                                  template_pdb = template_pdb,
                                                  ens_size = ens_size,
                                                  target_b = target_pdb.target_b[target_pdb.mask],
                                                  mask = target_pdb.mask,
                                                  rb_type = 'mix',
                                                  step = step,
                                                  start_step_list = start_step_list,
                                                  l = l)
                r = list(rb_optimizer.result)
                trans = r[0]
                rot = r[1:4]
                com = r[4:7]
            
    # check for errors
    else:
        raise Exception('No or wrong aniso input')

###########################################################################
#                       Analysis of ensemble                              #
###########################################################################

    l.process_message('Analysing final rigid body ensemble...')

    #### Final ensemble ####

    #template_pdb.er_ens_hierarchy
    final_b = template_pdb.rb_ens_B.as_numpy_array()
    l.show_info('Average CA B-factor: {:3.2f}'.format(final_b.mean()))
    l.show_info('Average rmsf (A)   : {:3.2f}'.format(b_to_rmsf(final_b.mean())))

    #### Plot target and fitted B-factors #####

    plot_b_factors(residue_numbers,target_pdb.target_b,final_b,p.output.plot_out)

    #### Plot histogram of translation vector lengths ####

    if p.input.rb_type != 'rot':
        # histogram of x translation
        plot_histo(np.array(template_pdb.v_list)[:,0],
                   'x_rb_trans_vectors.eps',
                   'X coord ($\AA$)')

        # histogram of y translation
        plot_histo(np.array(template_pdb.v_list)[:,1],
                   'y_rb_trans_vectors.eps',
                   'Y coord ($\AA$)')

        # histogram of z translation
        plot_histo(np.array(template_pdb.v_list)[:,2],
                   'z_rb_trans_vectors.eps',
                   'Z coord ($\AA$)')

        # histogram of magnitude of translation vectors
        plot_histo(np.array(template_pdb.v_dist_list),
                   'rb_trans_vectors.eps',
                   'Translation Vector Length ($\AA$)')

    #### Analysis of Rotation angles ####

    if p.input.rb_type != 'trans' and p.params.aniso == False:
        l.show_info('Average rotation angle (deg): {}'.format(np.mean(template_pdb.angle_list)))
        plot_histo(np.rad2deg(template_pdb.angle_list),
                   'rb_rot_angles.eps',
                   'angle ($^\circ$)')

###########################################################################
#                        Output of ensemble                               #
###########################################################################

    #### Writing RB ensemble to file

    template_pdb.rb_hierarchy.write_pdb_file(file_name=p.output.pdb_out,
                                             open_append=False,
                                             crystal_symmetry=target_pdb.symmetry,
                                             write_scale_records=True,
                                             append_end=False,
                                             interleaved_conf=0,
                                             atom_hetatm=True)

    l.show_info('Score: {}'.format(template_pdb.rb_score))

    #### End of Program ####

    return l

def b_to_rmsf(b):
    '''
    Returns dispalcement (u) for certain B-value
    '''

    return np.sqrt(b / (8 * (np.pi**2) ))

def init_plot():
    '''
    Initialize plot and global parameters
    '''

    #### Prep plot ####

    fig,ax = plt.subplots(1,1,figsize=(7.5,6))

    #### Adjust global parameters ####

    params = {
        'axes.labelsize' : 16,
        'legend.fontsize' : 14,
        'xtick.labelsize' : 12,
        'ytick.labelsize' : 12,}

    matplotlib.rcParams.update(params)

    #### Return results ####

    return fig,ax

def plot_layout(ax):
    '''
    Adjust axes / spines for plts
    '''

    #### Layout ####

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    return ax

def plot_histo(l,nm,x):
    '''
    Plot histogram of list
    '''

    #### Init plot ####

    fig,ax = init_plot()

    #### Plot histogram ####

    ax.hist(l,
            bins = 15,
            normed = True,
            color = (0.2,0.6,0.2))


    #### Add additional info ####

    ax.set_xlabel(x)
    ax.set_ylabel('Frequency (normed)')

    #### Add Layout ####

    ax = plot_layout(ax)

    #### Save plot ####

    fig.savefig(nm,
                format = 'eps',
                transparent=True)

    #### Clear figure ####

    fig.clf()

def plot_b_factors(x,y1,y2,nm):
    '''
    Simple plotting method

    - Plots 2 lines against same x-value list
    - Does some layouting
    - Writes to file (nm, should be .eps)
    '''

    fig,ax = init_plot()

    #### Plot lines ####

    ax.plot(x,y1,
            label = 'Target B-factor',
            lw = 3,
            c = (0,0,0))
    ax.plot(x,y2,
            label = 'Fitted B-factor',
            lw=2)

    #### Add additional info ####

    ax.legend()
    ax.set_xlabel('Residue')
    ax.set_ylabel('B-factor ($\AA^{2})$')
    ax.set_ylim((0,np.max(y1)))

    #### layout plot ####

    ax = plot_layout(ax)
    ax.grid()

    #### Write to file ####

    fig.savefig(nm,
                format = 'eps',
                transparent=True)

    #### Clear figure ####

    fig.clf()
