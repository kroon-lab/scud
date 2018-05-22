 #!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

from scitbx.array_family import flex

import phil as addLoes_phil
from scud.general.log_writer import Log
from scud.pdb.PDB import PDBClass
from scud.exp.EXP import EXPClass
from scud.mtz.slice_mtz import slice_mtz
from scud.plt.Plot import Plot
from scud.pdb.ens2b import ens2b

def run(args=None, l= None):
    '''
    jan, 2018
    Additivity check

    TODO:
    - 2 PDF files:
      - Single structre
      - Ensemble, internal motion only
    - Input sigma for translation is similar
    - Create 2 disorder models from same translation sigma's:
      - Translation only
      - Internal + Translation
    - Plot B-factor from ensembles:
      - Internal only
      - Internal + Translation
      - Translation only
    - Calculate diffuse scattering for 3 models
      - Internal only
      - Translation only
      - Internal + Translation
    - Slice IDFF 0kl,h0l,hk0 from mtz, save as np.array
    - Plots
    - Subtract Translation only array from Internal+Translation
    - Plot difference array
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='addLoes_log.txt')
    l.title("exp.addLoes module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = addLoes_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.addLoes

###########################################################################
#                           Input files                                   #
###########################################################################

    # If translation models should be created:
    if p.params.do_rb == True:
        l.process_message('Applying rigid body operation')
        # Read files
        single_pdb = PDBClass(fname=p.input.single_pdb, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        internal_pdb = PDBClass(fname=p.input.internal_pdb, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        l.process_message('Input models read...')
        # Input params
        sigma = p.params.trans_sigma
        nr_out = p.params.nr_models
        l.show_info('Tanslation sigma: {}'.format(sigma))

###########################################################################
#                           Create Models                                 #
###########################################################################
        
        # Translate single_pdb:
        single_pdb.model_list = single_pdb.create_random_model_list(l = l,
                                                                    ensemble_size = nr_out)
        single_pdb.set_B_zero()
        # Set occupancy to 1
        single_pdb.set_occ()
        center_of_mass =  single_pdb.calculate_center_of_mass()
        # Create final ensemble, ER result + added btls motion
        single_pdb.make_rot_trans_ens(trans_only = True,
                                    rot_only = False,
                                    trans_sigma = sigma,
                                    rot_sigma = 0.0,
                                    center_of_mass = center_of_mass,
                                    l = l)
        # Write final ensemble
        single_pdb.write_hierarchy_to_pdb(out_name=p.input.translation_pdb)
        l.process_message('Translation only done...')
        l.show_info('Written to: {}'.format(p.input.translation_pdb))
        
        # Translate internal internal_pdb
        internal_pdb.model_list = internal_pdb.create_random_model_list(l = l,
                                                                    ensemble_size = nr_out)
        internal_pdb.set_B_zero()
        # Set occupancy to 1
        internal_pdb.set_occ()
        center_of_mass =  internal_pdb.calculate_center_of_mass()
        # Create final ensemble, ER result + added btls motion
        internal_pdb.make_rot_trans_ens(trans_only = True,
                                    rot_only = False,
                                    trans_sigma = sigma,
                                    rot_sigma = 0.0,
                                    center_of_mass = center_of_mass,
                                    l = l)
        # Write final ensemble
        internal_pdb.write_hierarchy_to_pdb(out_name=p.input.int_translation_pdb)
        l.process_message('Translation + Internal motion done...')
        l.show_info('Written to: {}'.format(p.input.int_translation_pdb))

###########################################################################
#                           Plot B-factors                                #
###########################################################################
        
    # If true, run ens2b on the created ensembles
    if p.params.do_bfactor == True:
        fls = [p.input.translation_pdb,
               p.input.internal_pdb,
               p.input.int_translation_pdb]
        for f in fls:
            args = ['pdb_in={}'.format(f),
                    'plot_dat={}'.format(f[:-4]+'.json')]
            ens2b.run(args=args,
                      l = l)
            
###########################################################################
#                           Calculate IDFF                                #
###########################################################################

    # If diffuse scattering should be calculated from models
    if p.params.do_diff_calc == True:
        l.process_message('Performing diffuse scattering calculations...')
        # Parameters for diffuse scattering calculation
        supercell_num = 100
        size_h = p.params.size_h
        size_k = p.params.size_k
        size_l = p.params.size_l
        Ncpu = 10

        # Loop over three models
        for pdb in [p.input.internal_pdb,
                    p.input.translation_pdb,
                    p.input.int_translation_pdb]:
            l.show_info('PDB in: {}'.format(pdb))
            # Easy method to calculate diffuse scattering
            ex = EXPClass()    
            ex.calc_diffuse_and_map(pdb = pdb,
                                    supercell_num = supercell_num,
                                    size_h = size_h,
                                    size_k = size_k,
                                    size_l = size_l,
                                    Ncpu = Ncpu,
                                    l = l)
        l.process_message('Diffuse scattering calculated...')

###########################################################################
#                           Slice MTZ files                               #
###########################################################################

    # If true, slice IDFF mtz files, save 3 npy arrays
    if p.params.do_slice == True:
        l.process_message('Slicing mtz files and writing .npy arrays')
        # Input parameters
        fls = [p.input.translation_pdb[:-4]+'.mtz',
               p.input.int_translation_pdb[:-4]+'.mtz',
               p.input.internal_pdb[:-4]+'.mtz']
        # Plotting parameters
        vmin = p.params.vmin
        vmax = p.params.vmax

        # Loop over mtz files
        for i,n in enumerate(fls):
            if i == 2:
                vmin = p.params.vmin_internal
                vmax = p.params.vmax_internal
            l.show_info('mtz: {}'.format(n))
            args = ['mtz_in={}'.format(n),
                    'vmin={}'.format(vmin),
                    'vmax={}'.format(vmax),
                    'array=IDFF',
                    'write_slices=True']
            slice_mtz.run(args = args,
                          l = l)

###########################################################################
#                           Subtractions                                  #
###########################################################################

    import matplotlib.pyplot as plt
    # Perform actual experiment
    if p.params.do_subtractions == True:
        l.process_message('Subtracting Arrays')
        vmin = p.params.vmin_sub
        vmax = p.params.vmax_sub
        for ind in ['h','k','l']:
            total = np.load(p.input.int_translation_pdb[:-4]+'_{}.npy'.format(ind))
            s_1 = np.load(p.input.translation_pdb[:-4]+'_{}.npy'.format(ind))
            sub = total - s_1
            thing = ((sub > 0.0 ).astype(float) * 2) -1
            print thing[0,0]
            fig = plt.figure(figsize=(9,9))
            fig.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0,wspace=0.0,hspace=0.0)
            ax1 = fig.add_subplot(111)
            # Turn axes (x,y) off
            ax1.set_axis_off()
            # Start plot
            img = ax1.imshow(np.rot90(thing), vmin=-1, vmax=1)
            img.set_cmap('binary_r')
            # Force square plot
            ax1.set_aspect('equal')
            # Add colorbar
            colorbar = False
            fig.savefig('sum_binary_'+ind+'.eps',format='eps',transparent=True, bbox_inches='tight')
            for ar in [total,s_1,sub]:
                print 'mean: {:12.2f}, min: {:12.2f}, max: {:12.2}'.format(np.nanmean(ar),
                                                                           np.nanmin(ar),
                                                                           np.nanmax(ar))

            trans = np.load(p.input.translation_pdb[:-4]+'_{}.npy'.format(ind))
            internal = np.load('internal'+'_{}.npy'.format(ind))
            add = trans + internal
            pt = Plot()
            pt.contour2D(add,
                         plotfile_name='trans_int_sum_'+ind+'.eps',
                         vmin=1000, vmax=75000000)
            
