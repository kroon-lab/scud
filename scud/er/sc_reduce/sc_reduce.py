from __future__ import division
from cctbx import crystal
from cctbx import uctbx
from cctbx import miller
import iotbx.pdb
from iotbx import mtz
import mmtbx.utils
from libtbx import easy_mp
from scitbx.array_family import flex

import numpy as np
import os,sys
import gc

import phil as sc_reduce_phil
from scud.general.log_writer import Log

from scud.pdb.PDB import PDBClass
from scud.pdb.supercell import supercell

def run(args=None, l=None):
    '''
    * read in a supercell structure
    * fractionalize the structure based on the 'orginial' unit cell
    * Brute force find the best super cell translation and symmetry 
      operation per molecule to place it back in the asu
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='er.sc_reduce.txt')
    l.title("er.sc_reduce module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = sc_reduce_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.sc_reduce

    # Read input super cell PDB
    sc_pdb = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    # Read unit cell (super cell)
    sc_unit_cell = sc_pdb.symmetry.unit_cell()

    # Read reference PDB
    ref_pdb = PDBClass(fname=p.input.ref_pdb, selection="not (resname HOH) and not (resname CL) and not (resname NA)")
    ref_pdb.hierarchy.remove_alt_confs(always_keep_one_conformer=True)
    
    # Read reference unit cell
    ref_unit_cell = ref_pdb.symmetry.unit_cell()

###########################################################################
#             Find inverse symmetry operation and ref coordinates         #
###########################################################################

    l.process_message('Listing inverse symmetry operations...')

    #### Generating sym ops based on space group ####

    # List to save in
    inv_sym_op_ls = []
    
    # Loop over symmetry operation
    for sym_op in ref_pdb.symmetry.space_group():
        
        # Extract rotation and translation
        rot = sym_op.r().as_double()
        trans = sym_op.t().as_double()

        # convert to numpy
        rot = np.array(rot)
        trans = np.array(trans)

        # Print info
        l.show_info('t: {}'.format(trans))
        l.show_info('r: {}'.format(rot[0:3]))
        l.show_info('r: {}'.format(rot[3:6]))
        l.show_info('r: {}'.format(rot[6:9]))
        l.show_info('\n')

        # Invert symmetry operations
        inv_rot = tuple(np.transpose(np.array([rot[0:3],rot[3:6],rot[6:9]])).flatten())
        inv_trans = tuple(trans*-1)
        
        # Save to list
        inv_sym_op_ls.append([inv_rot,inv_trans])

    #### Extract reference coordinates ####

    # get model:
    ref_model = ref_pdb.detached_model_copy(n=0)
    
    # fractionalize coordinates:
    ref_model.atoms().set_xyz(ref_unit_cell.fractionalize(
        ref_model.atoms().extract_xyz()))

    #### Save reference coordinates (ONE CHAIN) ####
    
    # reference coordinates:
    ref_coords = ref_model.chains()[0].atoms().extract_xyz()
    
###########################################################################
#              Brute force reduce super cell to ASU                       #
###########################################################################

    #### Create final hierarchy ####

    final_hierarchy = iotbx.pdb.hierarchy.root()

    # Loop over all models in the super cell PDB
    for m,model in enumerate(sc_pdb.hierarchy.models()):
        
        # Detach model
        det_model = sc_pdb.detached_model_copy(n=m)
        
        # Fractionalize based on ref unit cell:
        det_model.atoms().set_xyz(
            ref_unit_cell.fractionalize(det_model.atoms().extract_xyz())
        )

        # Now per chain within model try all sc translations and symemtry operations
        for i,chain in enumerate(det_model.chains()):

            # Hacky to obtain an empty model object
            new_model = model.detached_copy()
            for c in range(len(new_model.chains())):
                new_model.remove_chain(0)

            # Possible super cell translations
            sc_trans_ls = sc_trans_list(p.input.supercell_size)

            #### Try all translations + inverse symmetry operations ####

            # Save coordinates belonging to (smallest) score
            final_coords = None
            score = None

            # Loop over all possible supercell translations
            for sc_trans in sc_trans_ls:

                # Loop over all symmetry ops
                for inv_sym_op in inv_sym_op_ls:

                    #### Generate possible answer ####
                    
                    # Deep copy chain
                    tst_chain = chain.detached_copy().atoms()

                    # Apply super cell translation
                    tst_chain.set_xyz(tst_chain.extract_xyz() + sc_trans)

                    # Apply inverse symmetry operation
                    tst_chain.set_xyz((tst_chain.extract_xyz()
                                       +inv_sym_op[1])*inv_sym_op[0]
                    )

                    #### Score possible solution ####
                    
                    # Calculate score
                    tst_score = calc_score(tst_chain,
                                           ref_coords)

                    # Get initial values
                    if score == None:
                        score = tst_score
                        final_coords = tst_chain.extract_xyz().deep_copy()
                    elif score != None:
                        # If sore is better:
                        if tst_score < score:
                            score = tst_score
                            final_coords = tst_chain.extract_xyz().deep_copy()
                        elif tst_score > score:
                            continue
                        elif tst_score == score:
                            continue
                        else:
                            raise Exception('error!!')
                    else:
                        raise Exception('Score is not correct type')

            l.show_info('model: {}, chain: {}'.format(m,i))

            #### Prep output ####
            
            # rename new chain (A because 1 chain per model) 
            chain.id = 'A'
            
            # Convert chain atoms back to cartesian
            chain.atoms().set_xyz(
                ref_unit_cell.orthogonalize(final_coords)
                )
            
            # Append chain to empty model
            new_model.append_chain(chain.detached_copy())
            # Append model to final hierarchy
            final_hierarchy.append_model(new_model)
            gc.collect(); gc.collect(); gc.collect()

###########################################################################
#                               Write to file                             #
###########################################################################
            
    # New space group / unit cell based on ref pdb (FOR NOW IT IS IN P1!!)
    new_sym=crystal.symmetry(ref_unit_cell.parameters(),
                             space_group= 'P1')

    # write to file:
    final_hierarchy.write_pdb_file(crystal_symmetry = new_sym,
                                   file_name = p.input.pdb_in[:-4]+'_asu.pdb')

    return l
                                

def calc_score(atom_list,
               ref_list):
    '''
    Calculate 'alignement score'

    all distances between same atoms summed

    '''
    return np.sum(flex.sqrt((atom_list.extract_xyz() - ref_list).dot()).as_numpy_array())


def sc_trans_list(i):
    '''
    Create a list of all possile combinations between -i and i+1 in (x,y,z)
    '''
    return np.array([(x,y,z)
                     for x in range(-1*i,i+1)
                     for y in range(-1*i,i+1)
                     for z in range(-1*i,i+1)])
    
                
# def old():
#     cm_sym_ls = []
#     for sym_op in ref_pdb.symmetry.space_group():
#         rot = sym_op.r().as_double()
#         trans = sym_op.t().as_double()
#         print rot, trans
#         det_model = ref_pdb.detached_model_copy(n=0)
#         # Fractionalize ref PDB
#         det_model.atoms().set_xyz(ref_unit_cell.fractionalize(
#             det_model.atoms().extract_xyz()))
#         # Place ASU within unit cell boundaries
#         det_model.atoms().set_xyz(det_model.atoms().extract_xyz()+check_cm(det_model.atoms().extract_xyz().mean()))
#         det_model.atoms().set_xyz(rot*det_model.atoms().extract_xyz()
#                                   +trans)
#         cm = det_model.atoms().extract_xyz().mean()
#         print cm
#         # det_model.atoms().set_xyz(det_model.atoms().extract_xyz()+check_cm(cm))
#         rot = np.array(rot)
#         inv_rot = tuple(np.transpose(np.array([rot[0:3],rot[3:6],rot[6:9]])).flatten())
#         inv_trans = tuple(np.array(trans)*-1)
#         cm_sym_ls.append([tuple(np.array(cm)+check_cm(cm)),inv_rot,inv_trans])

    
# ###########################################################################
# #                          Reduce to P1 and inverse symmetry              #
# ###########################################################################


#     l.process_message('Placing chains in P1 and then apply sym op to create ASU...')
#     # Final hierarchy, one model per chain
#     final_hierarchy = iotbx.pdb.hierarchy.root()
# #    help(flex.vec3_double)
    
#     xyz = flex.vec3_double([(1,0,1),(1,0,-1)])
    
#     for sym_op in ref_pdb.symmetry.space_group():
#         rot = sym_op.r().as_double()
#         trans = sym_op.t().as_double()
#         print rot, trans

#         print 'start', list(xyz)[0]
        
#         r = xyz*rot
        
#         print 'rot', list(r)[0]

#         t = r+trans
        
#         print 'trans', list(t)[0]
#         print '-----'

        
#     for i in cm_sym_ls: print i
#     for m,model in enumerate(sc_pdb.hierarchy.models()):
#         # Detach model
#         det_model = sc_pdb.detached_model_copy(n=m)
#         # Fractionalize based on ref unit cell:
#         det_model.atoms().set_xyz(
#             ref_unit_cell.fractionalize(det_model.atoms().extract_xyz())
#         )
        
#         # Now per chain within model find sym op and reduce to ASU
#         for chain in det_model.chains():
#             # Long winded hack to obtain an empty model object
#             new_model = model.detached_copy()
#             for c in range(len(new_model.chains())):
#                 new_model.remove_chain(0)

#             print type(chain.atoms().extract_xyz())
                
#             #### CoM of the chain ####
            
#             cm_chain = chain.atoms().extract_xyz().mean()

#             print '->', cm_chain

#             #### Find best CoM correction and the symmetry index closest to this ####
            
#             cm_correction, sym_index = dist_calc(cm_chain,cm_sym_ls)


#             chain.atoms().set_xyz(chain.atoms().extract_xyz()-cm_correction)

#             print 'sc_trans', cm_correction

#             print 'inverse unit cell trans', cm_sym_ls[sym_index][2]

#             print 'inverse unit cell rot', cm_sym_ls[sym_index][1]
            
#             chain.atoms().set_xyz((chain.atoms().extract_xyz()
#                                    +cm_sym_ls[sym_index][2])*cm_sym_ls[sym_index][1]
#             )

#             new_cm = chain.atoms().extract_xyz().mean()

#             print 'new cm', new_cm

#             if new_cm[0] < 0.0 or new_cm[0] > 1.0 or new_cm[1] < 0.0 or new_cm[1] > 1.0 or new_cm[2] < 0.0 or new_cm[2] > 1.0:
#                 continue

            
#             # correct_bool, sym_op_index = inv_op_check(cm_uncorrected,
#             #                                           cm_corrected,
#             #                                           cm_sym_ls)

#             # print correct_bool, sym_op_index

#             # correct_bool = True
            
#             # if correct_bool == True:
#             #     chain.atoms().set_xyz(chain.atoms().extract_xyz()+check_cm(cm_uncorrected))

#             # Apply inverse symmetry operation
#             # chain.atoms().set_xyz((chain.atoms().extract_xyz()
#             #                        +cm_sym_ls[sym_op_index][2])*cm_sym_ls[sym_op_index][1]
#             # )

#             # quit()
            
#             # # Check cm per chain and correct if not in P1 cell
#             # cm = chain.atoms().extract_xyz().mean()
#             # op_num = inv_op(cm,
#             #                 cm_sym_ls,
#             #                 check = False)
#             # chain.atoms().set_xyz(chain.atoms().extract_xyz()+check_cm(cm))
#             # # Find sym op based on ref locations and correct
#             # new_cm = chain.atoms().extract_xyz().mean()
#             # op_num = inv_op(new_cm,
#             #                 cm_sym_ls,
#             #                 check=True)
#             # # Apply inverse symmetry operation
#             # chain.atoms().set_xyz((chain.atoms().extract_xyz()
#             #                        +cm_sym_ls[op_num][2])*cm_sym_ls[op_num][1]
#             # )
            
#             # Double check if everything is within the unit cell
#             # cm = chain.atoms().extract_xyz().mean()
#             # chain.atoms().set_xyz(chain.atoms().extract_xyz()+check_cm(cm))
#             # Rename chain
#             chain.id = 'A'
#             # Convert chain atoms back to cartesian
#             chain.atoms().set_xyz(
#                 ref_unit_cell.orthogonalize(chain.atoms().extract_xyz())
#                 )
#             # Append chain to empty model
#             new_model.append_chain(chain.detached_copy())
#             # Append model to final hierarchy
#             final_hierarchy.append_model(new_model)
#     # New space group / unit cell based on ref pdb
#     new_sym=crystal.symmetry(ref_unit_cell.parameters(),
#                              space_group= 'P1')
#     # write to file:
#     final_hierarchy.write_pdb_file(crystal_symmetry = new_sym,
#                                    file_name = p.input.pdb_in[:-4]+'_asu.pdb')
    
#     return l

# def inv_op_check(cm_uncorrected,
#                  cm_corrected,
#                  cm_sym_ls):
#     '''
#     input:
#       - corrected CoM (placed within P1 fractionalized cell)
#       - uncorrected CoM (Just outside)
#       - Ref CoM+Inv_Symop list

#     Find if corrected or uncorrected-1 has the shortest distance to a reference molecule

#     return if there should be a CoM correction and the inverse_sym_op_list_index
#     '''

#     print cm_uncorrected, cm_corrected
    
#     d_1, sym_index_1 = dist_calc(cm_uncorrected,cm_sym_ls)

#     if (d_1 - 1) < d_2:
#         return False, sym_index_1
#     elif (d_1 - 1) > d_2:
#         return True, sym_index_2
#     elif (d_1 - 1) == d_2:
#         return True, sym_index_2
#     else:
#         raise Exception('Error in finding symmetry operation')
    
# def dist_calc(cm,ref):
#     '''
#     Try back translations compared to ref CoM to find which molecule is similar to which
#     '''

#     cm_corr_ls = np.array([[x,y,z]
#                            for x in range(-2,3)
#                            for y in range(-2,3)
#                            for z in range(-2,3)])

#     # Will contain d,cm_corr,inv_sym_index
#     d_ls = []
#     inv_sym_index_ls = []
#     corr_ls = []
#     for cm_corr in cm_corr_ls:
#         d_tmp = []
#         for c in ref:
#             tst_cm = cm - cm_corr
#             d = np.linalg.norm(tst_cm - np.array(c[0]))
#             d_tmp.append(d)
#         d_ls.append(np.sort(d_tmp)[0])
#         inv_sym_index_ls.append(np.argmin(d_tmp))
#         corr_ls.append(cm_corr)
#     ind = np.argmin(d_ls)
#     return tuple(corr_ls[ind]), inv_sym_index_ls[ind]

# def inv_op(cm,
#            cm_sym_ls,
#            check = False):
#     '''
#     Find the closest reference center of mass
#     '''
#     cm = np.array(cm)
#     diff_ls = []
#     for c in cm_sym_ls:
#         cc = np.array(c[0])
#         d = np.linalg.norm(cm - cc)
#         diff_ls.append(d)
#     if check == False:
#         print '===='
#     elif check == True:
#         print '++++'
#     print np.sort(diff_ls)[0:2]
#     cm_index = np.argmin(diff_ls)
#     return cm_index
    
# def check_cm(cm):
#     '''
#     Find translation vector for center of mass

#     This only works for 1x1x1 and 2x2x2 cells!!!
#     '''

#     # (2.2, -1.1, 0.1) -> (-2, -1, 0)
#     t = []
#     for i in cm:
#         if 0.0 < i < 1.0:
#             t.append(0)
#         if i < 0.0:
#             if -1.0 < i < 0.0:
#                 t.append(1)
#             else:
#                 t.append((abs(i)//1))
#         if i > 1.0:
#             t.append((i//1)*-1.)
            
#     # n_cm = []
    
#     # for i in cm:
#     #     if 0.0 < i < 1.0:
#     #         t.append(0)
#     #     if i < 0.0:
#     #         t.append(
#     # t = []
#     # for i in cm:
#     #     if 0.0 < i < 1.0:
#     #         t.append(0)
#     #     if i < 0.0:
#     #         t.append(1)
#     #     if i > 1.0:
#     #         t.append(-1)
#     return tuple(t)
