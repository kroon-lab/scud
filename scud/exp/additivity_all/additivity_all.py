#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os, glob, shutil

import libtbx.phil
import phil as additivity_all_phil

from scud.general.log_writer import Log
from scud.general.phil_methods import init_command_line_phil
from scud.general.method_library import change_element_in_json

from scud.pdb.rigid_body import rigid_body
from scud.pdb.supercell import supercell
from scud.pdb.ens2b import ens2b

from scud.mtz.slice_mtz import slice_mtz
from scud.mtz.mtz2map import mtz2map

from scud.exp.additivity import additivity

from scud.plt.lines import lines

def run(args=None):
    '''
    Check additivity of motion in reciprocal space, subtract mtz files
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='additivity_all_log.txt')
    l.title("exp.additivity_all module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = additivity_all_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.additivity_all
    
###########################################################################
#                        Generate ensembles                               #
###########################################################################

    ens_list = ['rbTrans','rbRot','rbMix','ens','rbTransEns','rbRotEns','rbMixEns']
    ens_command_list = [[True,['pdb_in={}'.format(p.input.single_pdb),
                               'pdb_out={}.pdb'.format(ens_list[0]),
                               'translation_only=True']],
                        [True,['pdb_in={}'.format(p.input.single_pdb),
                               'pdb_out={}.pdb'.format(ens_list[1]),
                               'rotation_only=True']],
                        [True,['pdb_in={}'.format(p.input.single_pdb),
                               'pdb_out={}.pdb'.format(ens_list[2])]],
                        [False,[p.input.ensemble_pdb]], # Only ensemble no rbMotion
                        [True,['pdb_in={}'.format(p.input.ensemble_pdb),
                               'pdb_out={}.pdb'.format(ens_list[4]),
                               'translation_only=True']],
                        [True,['pdb_in={}'.format(p.input.ensemble_pdb),
                               'pdb_out={}.pdb'.format(ens_list[5]),
                               'rotation_only=True']],
                        [True,['pdb_in={}'.format(p.input.ensemble_pdb),
                               'pdb_out={}.pdb'.format(ens_list[6])]]]

    if p.params.calc_diff:
        batch_diffuse_calc(p = p,
                           ens_list = ens_list,
                           ens_command_list=ens_command_list,
                           l = l)
###########################################################################
#                       Create big slice image                            #
###########################################################################

    if p.params.merge_slices: merge_slices(ens_list = ens_list,
                                           l = l)

###########################################################################
#                         Plot all B-Factors                              #
###########################################################################

    if p.params.plot_all_B_factors: plot_all_B_factors(ens_list = ens_list,
                                                       l = l)

###########################################################################
#                   Subtract all scattering stuff                         #
###########################################################################

    # All sums, per element [0] - [1] = [2]
    sum_list = [['rbMix','rbRot','rbTrans'],
                ['rbMix','rbTrans','rbRot'],
                ['rbTransEns','rbTrans','ens'],
                ['rbTransEns','ens','rbTrans'],
                ['rbRotEns','rbRot','ens'],
                ['rbRotEns','ens','rbRot'],
                ['rbMixEns','rbMix','ens'],
                ['rbMixEns','ens','rbMix']]
    # special sum: [0] - [1] - [2] = [3]
    special = ['rbMixEns','rbTrans','rbRot','ens']

    # List to be appended to ens_list for all R and CC value calculations
    outname_list = []
    # Loop over all summations
    for sum in sum_list:
        # File name
        outname = '{}_{}'.format(sum[0],sum[1])
        outname_list.append(outname)
        # Subtract sum[1] from sum[0]
        a = ['mtz_1={}_diffuse.mtz'.format(sum[0]),
             'mtz_2={}_diffuse.mtz'.format(sum[1]),
             'mtz_out={}_diffuse.mtz'.format(outname)]
        additivity.run(args = a,
                       l = l)
        # Create 2D slice and map from new mtz
        slice_and_make_map(nm = outname,
                           l = l)
        # Combine slices to final png file
        os.system('convert {}_slice.eps {}_slice.eps {}_slice.eps {}_slice.eps  +append {}_slice.png'.format(sum[0], sum[1],outname,sum[2],outname))
        
###########################################################################
#                   Calculate all R and CC values                         #
###########################################################################

    # use mtz.correlation to calculate values
    l.close_log()


def slice_and_make_map(nm = None,
                       l = None):
    '''
    Slice a mtz to make 2D plot, and create a 3D map using
    slice_mtz and mtz2map respectively.
    '''
    a = ['mtz_in={}_diffuse.mtz'.format(nm),
         'array=IDFF',
         'mtz_out={}_slice.mtz'.format(nm),
         'npy_out={}_slice.npy'.format(nm),
         'plt_out={}_slice.eps'.format(nm)]
    slice_mtz.run(args=a,
                      l=l)
    # Calculate map from mtz
    a = ['mtz_in={}_diffuse.mtz'.format(nm),
         'array=IDFF',
         'map_out={}.map'.format(nm)]
    mtz2map.run(args=a,
                l = l)
    
def plot_all_B_factors(ens_list = None,
                       l = None):
    '''
    use plt.lines to plot the B-factors calculated from the ensemble rmsf's
    '''
    # having all json file with ensemble B-factors use plt.lines
    json_dat_list = ['dat={}_ensembleB.json'.format(i) for i in ens_list]
    json_dat_list.append('plot_out=all_ens_B_lines.eps')
    a = json_dat_list
    lines.run(args = a,
              l = l)
    return l
    
def merge_slices(ens_list = None,
                 name = '_slice',
                 l=None):
    '''
    Use 'convert' a gimp plugin to merge all slice images into 1 image
    '''
    # Merge all slice eps images into one file 
    l.process_message('Merging slice eps images...')
    eps_string_1 = ''
    eps_string_2 = ''
    n = int(len(ens_list) / 2)
    for s in ens_list[:n]: eps_string_1 += ' {}{}.eps'.format(s,name)
    for s in ens_list[n:]: eps_string_2 += ' {}{}.eps'.format(s,name)
    os.system('convert {} +append tmp_1.png'.format(eps_string_1))
    os.system('convert {} +append tmp_2.png'.format(eps_string_2))
    os.system('convert tmp_1.png tmp_2.png -append all_slices.png')
    os.system('convert all_slices.png scale 10 smaller_all_slices.png')

    
def batch_diffuse_calc(p = None,
                       ens_list = None,
                       ens_command_list=None,
                       l = None):
    '''
    For a list of ensemble names and their rbCommand list, generate ensembles
    claculate diffuse scattering, organize files properly
    '''
    l.process_message('Generating ensembles, calculating diffuse scattering \n Plotting slice and calculating diffuse maps for {} ensembles...'.format(len(ens_list)))
    for i,ens_name in enumerate(ens_list):
        # Generate input arguments
        l.process_message('Creating ensemble: {}'.format(ens_name))
        if ens_command_list[i][0] == True:
            a = ens_command_list[i][1]
            # Generate the ensemble
            rigid_body.run(args=a,
                           l=l)
        else:
            shutil.copyfile(p.input.ensemble_pdb,ens_name+'.pdb')
            a = ['pdb_in={}.pdb'.format(ens_name),'plot_dat=plot_ensemble_B.json']
            ens2b.run(args=a,l=l)
        # rename B-factor Json file
        json_file_name = '{}_ensembleB.json'.format(ens_name)
        os.rename('plot_ensemble_B.json',json_file_name)
        # change label in B-factor Json file
        change_element_in_json(json_file = json_file_name,
                               label_name = 'label',
                               new_element = '{}'.format(ens_name),
                               l = l)

###########################################################################
#                          Calculate diffuse                              #
###########################################################################

        l.process_message('Calculating diffuse scattering using supercell method...')
        a = ['pdb_in={}.pdb'.format(ens_name),
             'supercell_num=100']
        supercell.run(args=a,
                      l=l)
        # Rename output files and remove temporary files
        os.rename('result.mtz','{}_diffuse.mtz'.format(ens_name))
        for f in glob.glob('tmp_*.mtz'): os.remove(f)

###########################################################################
#                   Plot slice and calculate map                          #
###########################################################################

        # Slice 2D to generate plot
        l.process_message('Calculating 2D slice and 3D map...')
        slice_and_make_map(nm = ens_name,
                       l = l)    
