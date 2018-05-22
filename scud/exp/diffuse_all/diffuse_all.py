#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os,glob,shutil

# General CCTBX modules
import libtbx.phil
import phil as diffuse_all_phil
# Log writer
from scud.general.log_writer import Log
# Scud modules
from scud.pdb.PDB import PDBClass
from scud.pdb.rigid_body import rigid_body
from scud.pdb.b import b
from scud.pdb.ens2b import ens2b
from scud.plt.lines import lines
from scud.map.MAP import MAPClass
from scud.map.MAP import MapMapClass
from scud.map.map2map import map2map
# Experiment function (for doing a lot of calculations)
from scud.exp.EXP import EXPClass

def run(args=None):
    '''
    - Create a large amount of dynamic models, 
    - Calculate diffuse scattering using supercells
    - Convert IDFF column to maps
    - In depth analysis between all maps and (possibly) data
    - Return CC-table
    - Return difference projection table
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='diffuse_all_log.txt')
    l.title("scud.exp.diffuse_all module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = diffuse_all_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.diffuse_all
    
###########################################################################
#                           Create models                                 #
###########################################################################

    model_names = ['trans.pdb',
                   'rot.pdb',
                   'mix.pdb',
                   'int.pdb',
                   'intTrans.pdb',
                   'intRot.pdb',
                   'intMix.pdb']
    # PDB stuff
    single = p.input.single_pdb
    l.process_message('Cleaning up ensemble PDB...')
    ensemble = 'int.pdb'
    # clean up ensemble PDB
    ens_tmp = PDBClass(fname=p.input.ensemble_pdb, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    ens_tmp.set_occ()
    ens_tmp.set_B_zero()
    ens_tmp.write_hierarchy_to_pdb(output_hierarchy = ens_tmp.hierarchy,
                                   out_name = ensemble)
    if p.params.create_models:
        # Use scud.pdb.rigid_body to create models from single structure
        # and ensemble
        # Loat single PDB to find average B-factor needed for model generation
        l.process_message('Creating dynamic models (pdbs), as input for diffuse scattering calculation...')
        single_pdb = PDBClass(fname=single, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
        average_B = single_pdb.average_B()
        # Rigid Body commands
        wilson_B= 20.18
        rigid_body_commands = [['pdb_in={}'.format(single),
                                'translation_only=True',
                                'pdb_out=trans.pdb',
                                'b_ext={}'.format(wilson_B)],
                               ['pdb_in={}'.format(single),
                                'rotation_only=True',
                                'pdb_out=rot.pdb',
                                'b_ref={}'.format(single)],
                               ['pdb_in={}'.format(single),
                                'pdb_out=mix.pdb',
                                'b_ref={}'.format(single)],
                               ['pdb_in={}'.format(ensemble),
                                'pdb_out=intTrans.pdb',
                                'b_ext={}'.format(wilson_B)],
                               ['pdb_in={}'.format(ensemble),
                                'rotation_only=True',
                                'pdb_out=intRot.pdb',
                                'b_ext={}'.format(wilson_B)],
                               ['pdb_in={}'.format(ensemble),
                                'pdb_out=intMix.pdb',
                                'b_ext={}'.format(wilson_B)]]
        l.process_message('Creating models...')
        for a in rigid_body_commands:
            rigid_body.run(args = a,
                           l = l)
        
###########################################################################
#                           B-Factor Plot                                 #
###########################################################################

    if p.params.B_factor_plots:
        l.process_message('Calculating and plotting B-Factors from all models...')
        # Create a B-Factor plot of all models and save all as a json file
        # Classic B-factor
        b_argument = ['pdb_in={}'.format(single),
                      'plot_dat=single.json']
        b.run(b_argument)
        # Ensemble B-factor
        ens2b_argument = ['pdb_in={}'.format(ensemble),
                          'plot_dat=int.json']
        ens2b.run(ens2b_argument)
        # New models B-Factor
        for n in model_names:
            ens2b_argument = ['pdb_in={}'.format(n),
                              'plot_dat={}.json'.format(n[:-4])]
            ens2b.run(ens2b_argument)
        plt_list = ['dat={}.json'.format(n[:-4]) for n in model_names]
        plt_list = ['dat=single.json',
                    'dat=int.json'] + plt_list + ['plot_out=all_b_factors.eps',
                                                  'show=False']
        lines.run(plt_list)
        l.process_message('All B-factors plotted in all_b_factors.eps...')
    
###########################################################################
#                   Calculate Diffuse scattering                          #
###########################################################################

    # Calculate ITOT,IBRG and IDFF for each model in supercell space and
    # convert to IDFF maps
    ex = EXPClass()
    if p.params.calc_diff:
        l.process_message('Calculating diffuse scattering for all models...')
        for pdb in model_names:
            ex.calc_diffuse_and_map(pdb = pdb,
                                    supercell_num = p.params.supercell_num,
                                    size_h = p.params.size_h,
                                    size_k = p.params.size_k,
                                    size_l = p.params.size_l,
                                    Ncpu = p.params.Ncpu,
                                    l = l)
        
###########################################################################
#                           Create CC table                               #
###########################################################################

    # If a data map is present append this to the map list
    # If not just calculate all possible CC values and save to a big file
    map_list = ['{}_IDFF.map'.format(i[:-4]) for i in model_names]
    if p.params.calc_CC:
        l.process_message('Calculating CC values between all maps and / or data...')
        if p.input.data_map != None:
            map_list.append(p.input.data_map)
            # List with all map names
        # ADD FOR LOOPS
        # List where CC values will be stored in
        CC_list = np.zeros((len(map_list),len(map_list)))
        CC_brg_list = np.zeros((len(map_list),len(map_list)))
        CC_aniso_list = np.zeros((len(map_list),len(map_list)))
        for i in range(len(map_list)): CC_list[i,i] = 1.0
        for i in range(len(map_list)): CC_brg_list[i,i] = 1.0
        for i in range(len(map_list)): CC_aniso_list[i,i] = 1.0
        print CC_list
        print CC_brg_list
        print CC_aniso_list
        # Create double loop, map2map for every map-map combination
        l.process_message('looping over all maps...')
        for i_map1, map1 in enumerate(map_list):
            for i_map2, map2 in enumerate(map_list):
                # Do not do double calculations
                if i_map1 <= i_map2: continue
                l.show_info('{} vs. {}'.format(map1,map2))
                # Check for experimental map
                m1_dat, m2_dat = False, False
                if map1 == p.input.data_map: m1_dat = True
                if map2 == p.input.data_map: m2_dat = True
                # Arguments for map2map
                map2map_args = ['map_1={}'.format(map1),
                                'map_2={}'.format(map2),
                                'map_1_exp={}'.format(m1_dat),
                                'map_2_exp={}'.format(m2_dat),
                                'size_h={}'.format(p.params.size_h),
                                'size_k={}'.format(p.params.size_k),
                                'size_l={}'.format(p.params.size_l),]
                l = map2map.run(args = map2map_args,
                                l = l)
                # ADD FOR LOOPS
                # Find back CC in map2map output and append to table
                CC_fname = 'map2map_{}_{}/CC.txt'.format(map1[:-4],map2[:-4])
                with open(CC_fname,'r') as f:
                    CC = float(f.readlines()[0])
                CC_list[i_map1,i_map2] = CC

                # Find back CC_bragg in map2map output and append to table
                CC_fname = 'map2map_{}_{}/CC_brg.txt'.format(map1[:-4],map2[:-4])
                with open(CC_fname,'r') as f:
                    CC = float(f.readlines()[0])
                CC_brg_list[i_map1,i_map2] = CC

                # Find back CC_inter_bragg in map2map output and append to table
                CC_fname = 'map2map_{}_{}/CC_aniso.txt'.format(map1[:-4],map2[:-4])
                with open(CC_fname,'r') as f:
                    CC = float(f.readlines()[0])
                CC_aniso_list[i_map1,i_map2] = CC

###########################################################################
#                      Writing CC's to table                              #
###########################################################################

        # ADD FOR LOOPS

        # Convert file names
        map_names = [i[:-4] for i in map_list]
        cc_lines = []
        # Make column header line
        column_header = [''] + map_names
        ch = ''
        for i in column_header: ch += '{:>15s}'.format(i)
        cc_lines.append(ch)
        # Make other lines
        for i,li in enumerate(CC_list):
            # Create first column (name)
            line = '{:15s}'.format(map_names[i])
            # add CC values to line
            for c in li:
                line += '{:15.2f}'.format(c)
            cc_lines.append(line)
        # Write to file
        with open('all_CC.txt','w') as f:
            for line in cc_lines:
                print >> f, line

        # Convert file names
        map_names = [i[:-4] for i in map_list]
        cc_lines = []
        # Make column header line
        column_header = [''] + map_names
        ch = ''
        for i in column_header: ch += '{:>15s}'.format(i)
        cc_lines.append(ch)
        # Make other lines        
        for i,li in enumerate(CC_brg_list):
            # Create first column (name)
            line = '{:15s}'.format(map_names[i])
            # add CC values to line
            for c in li:
                line += '{:15.2f}'.format(c)
            cc_lines.append(line)
        # Write to file
        with open('all_CC_brg.txt','w') as f:
            for line in cc_lines:
                print >> f, line

        # Convert file names
        map_names = [i[:-4] for i in map_list]
        cc_lines = []
        # Make column header line
        column_header = [''] + map_names
        ch = ''
        for i in column_header: ch += '{:>15s}'.format(i)
        cc_lines.append(ch)
        # Make other lines        
        for i,li in enumerate(CC_aniso_list):
            # Create first column (name)
            line = '{:15s}'.format(map_names[i])
            # add CC values to line
            for c in li:
                line += '{:15.2f}'.format(c)
            cc_lines.append(line)
        # Write to file
        with open('all_CC_aniso.txt','w') as f:
            for line in cc_lines:
                print >> f, line

###########################################################################
#                Make a difference projection table                       #
###########################################################################
    
    # All difference projections into one big 'image' sorted same as CC table
    l.process_message('Generating difference projection table...')
    # Create stings with image names per row
    s_list = []
    for i_map1, map1 in enumerate(map_list):
        s = ''
        for i_map2, map2 in enumerate(map_list):
            if i_map1 <= i_map2: continue
            projection_fname = 'map2map_{}_{}/projection_difference_map.eps'.format(map1[:-4],map2[:-4])
            s += ' ' + projection_fname
        s_list.append(s)

    # Create rows
    f_list = []
    for i,s in enumerate(s_list):
        f = '{}_tmp.png'.format(i)
        os.system('convert {} +append {}'.format(s,f))
        f_list.append(f)

    # Create final image
    f_string = ''
    for f in f_list: f_string += ' ' + f
    os.system('convert {} -append all_projection_differences.png'.format(f_string))
