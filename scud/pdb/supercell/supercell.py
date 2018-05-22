from __future__ import division
from cctbx import crystal
from cctbx import uctbx
from cctbx import miller
import iotbx.pdb
from iotbx import mtz
import mmtbx.utils
from libtbx import easy_mp
import gc

import numpy as np
from random import randint
import os,sys

import phil as supercell_phil
from scud.pdb.PDB import PDBClass
from scud.general.log_writer import Log

def calc_average(args):
    '''
    For mtz's in args list, calculate <I> and <F>, returns [<F>,<I>]
    '''
    print('First averaging step')
    mtz_names = args
    # Catch input having length of 1
    if len(mtz_names) == 1:
        mtz_object = mtz.object(args[0])
        miller_arrays = mtz_object.as_miller_arrays()
        return [miller_arays[0],miller_arrays[0].as_intensity_array().data()]
    # Loop over given mtz files, first file serves as template, sum
    # f and i
    else:
        for c,mtz_name in enumerate(mtz_names):
            if c == 0:
                template_object = mtz.object(mtz_name)
                template_miller_arrays = template_object.as_miller_arrays()
                f = template_miller_arrays[0]
                i = f.as_intensity_array().data()
            else:
                mtz_object = mtz.object(mtz_name)
                f_temp = mtz_object.as_miller_arrays()[0]
                i_temp = f_temp.as_intensity_array().data()
                f += f_temp
                i += i_temp
                # Average f and i and return
                # Delete objects to save memory?
                del f_temp
                del i_temp
                gc.collect()
        f_average = f / len(args)
        i_average = i / len(args)
        gc.collect(); gc.collect(); gc.collect()
        return [f_average,i_average]

def diff_calc(args):
    '''
    Average final F's and I's, calculate Itot and Ibrg, and use 
    this to calculate Idiffuse, and write this to a file.
    '''
    print('Final Averaging step, and write final mtz')
    # Loop over i_av's and f_av's, use first as templates
    for c,fi in enumerate(args):
        if c == 0:
            f = fi[0]
            i = fi[1]
        else:
            f+=fi[0]
            i+=fi[1]
    # Final averaging step
    f_average = f / len(args)
    i_average = i / len(args)
    # Convert f_average to I_brg
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

def apply_operation(model=None, operation=None, unit_cell=None):
    '''
    Apply symmetry operation (unit cell then supercell)  on detached model
    '''
    rotation = operation[2].r().as_double()
    translation = operation[2].t().as_double()
    sc_translation = operation[1]
    # Normal Symmetry operation within unit cell then supercell translation
    model.atoms().set_xyz(
        ((rotation*model.atoms().extract_xyz())    # Apply rotation (cartesian)
        + unit_cell.orthogonalize(translation))    # Apply translation part of sym op (cartesian)
        + unit_cell.orthogonalize(sc_translation)) # Apply sc_translation (cartesian)
    return model

def supercell(arg):
    '''
    Multiprocess function

    Construct list containing all symmetry operations to fill a supercell
    of a certain size. Construct is, FFT it and write out structure factors.
    '''
    p = arg[0]
    l = arg[1]
    outfile_num = arg[2]
    print 'supercell number: {}'.format(outfile_num)
    pdb_out = p.output.pdb_out[:-4]+'_{}'.format(outfile_num)+'.pdb'
    mtz_out = 'tmp_{}'.format(outfile_num)+'.mtz'
    
###########################################################################
#                  Read and fix input ensemble                            #
###########################################################################

    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    pdb_object.renum_and_count_models()
    unit_cell = pdb_object.symmetry.unit_cell()
    sym_list = []
    for sym_op in pdb_object.symmetry.space_group(): sym_list.append(sym_op)
    sym_list = tuple(sym_list)

###########################################################################
#                  Generate translation/rotation lists                    #
###########################################################################

    range_h = range(0,p.input.size_h)
    range_k = range(0,p.input.size_k)
    range_l = range(0,p.input.size_l)
    # Generate translation list:
    translation_list =  [[h,k,l] for h in range_h
                         for k in range_k
                         for l in range_l]
    # Total number of operation
    operation_num = len(sym_list) * len(translation_list)
    # Randomly pick a model num the amount of operations
    model_assign_list = [randint(0,pdb_object.model_num) for i in xrange(0,operation_num)]
    # put everything together
    all_operations_list,m = [],0
    for trans in translation_list:
        for sym in sym_list:
            all_operations_list.append([model_assign_list[m],trans,sym])
            m+=1
    all_operations_list = sorted(all_operations_list)

###########################################################################
#                  Apply all symmetry operations                          #
###########################################################################

    pdb_object.set_B_zero()
    pdb_object.start_or_attach_hierarchy()
    for num,operation in enumerate(all_operations_list):
        detached_model = pdb_object.detached_model_copy(n=operation[0])
        new_model = apply_operation(model=detached_model,
                                    operation=operation,
                                    unit_cell = unit_cell)
        if num == 0:
            pdb_object.start_or_attach_hierarchy(detached_model=new_model)
        else:
            pdb_object.new_hierarchy.models()[0].transfer_chains_from_other(new_model)
    pdb_object.rechain(model=pdb_object.new_hierarchy.models()[0])

###########################################################################
#                  Finishing up, write PDB                                #
###########################################################################

    # Generating new unit cell
    uc_param = unit_cell.parameters()[0:3]
    new_cell = np.array([uc_param[0]*p.input.size_h,
                         uc_param[1]*p.input.size_k,
                         uc_param[2]*p.input.size_l])
    new_symmetry = crystal.symmetry(
        unit_cell = tuple(np.append(new_cell,
                                    np.array(unit_cell.parameters())[3:7])),
        space_group_symbol="P1")
    # outputting structures, format depends on input
    if p.params.write_pdb:
        print('Writing supercell to PDB file...')
        pdb_object.write_hierarchy_to_pdb(symmetry=new_symmetry,
                                          out_name=pdb_out)
    if p.params.write_cif:
        print('Writing supercell to CIF file...')
        pdb_object.write_hierarchy_to_cif(symmetry=new_symmetry,
                                          out_name=pdb_out.replace('.pdb','.cif'))

###########################################################################
#                 FT supercell hierarchy and write mtz                    #
###########################################################################

    if p.params.FFT_structure:
    # convert hierarchy to Xray structure
        xray_structure = pdb_object.new_hierarchy.extract_xray_structure(
            crystal_symmetry = new_symmetry)
    # Get default fmodel parameters
        params = mmtbx.command_line.fmodel.\
                 fmodel_from_xray_structure_master_params.extract()
    # Input high resolution
        params.high_resolution = p.params.high_resolution
    # Convert structure to structurefactors
        mmtbx.utils.fmodel_from_xray_structure(
            xray_structure = xray_structure,
            f_obs          = None,
            add_sigmas     = False,
            params         = params,
            twin_law       = None,
            twin_fraction  = None,
            out            = sys.stdout).write_to_file(file_name = mtz_out,
                                                       obs_type=complex)
        gc.collect(); gc.collect(); gc.collect()
        return mtz_out

def run(args=None, l=None):
    '''
    Generate supercells, using an ensemble in PDB format as input.
    Will find the symmetry operations belonging to the space group, and 
    the translatinos needed to fill a supercell. 
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='supercell_log.txt')
    l.title("pdb.supercell module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = supercell_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.supercell
    ss_num=p.input.supercell_num

###########################################################################
#                       Setup parallel processing                         #
###########################################################################

    # Create argument list for mp pool

    args = [[p,l,i] for i in range(ss_num)]
    # Start multiprocessor calculation of supercells and structurefactors
    result = easy_mp.pool_map(
        func=supercell,
        args=args,
        processes=p.params.Ncpu)

###########################################################################
#                   Split SF's depending on the amount                    #
###########################################################################

    # Average MTZ files in a parallel way:
    if len(result) <= 10:
        chunk_size = len(result)
    else:
        chunk_size = int(len(result)/10)
    # Setup argument list for mp averaging
    chunk_list = [result[i:i+chunk_size] for i in range(0,ss_num,chunk_size)]
    args = chunk_list
    # Start mp averaging
    result = easy_mp.pool_map(
        func = calc_average,
        args=args)
    # Final averaging step
    diff_calc(result)
    
    return l

