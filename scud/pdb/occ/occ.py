from __future__ import division
import iotbx.pdb

import numpy as np
from random import randint
import os,sys

import phil as occ_phil
from scud.pdb.PDB import PDBClass
from scud.general.log_writer import Log

def run(args):
    '''
    Read ensemble (PDB format) calculate average structure and output B-Factor
    calculated from the RMSF
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    l = Log(log_name='occ_log.txt')
    l.title('occ module')

###########################################################################
#                      Process input Params                               #
###########################################################################

# use phil to process input
    working_params = occ_phil.phil_parse(args=args,log=l)
    p = working_params.occ

###########################################################################
#                          Read  pdb_in                                   #
###########################################################################

    l.process_message('Reading and analyzing PDB file / Structure...')
    sel = p.params.occ_sel
    pdb_object = PDBClass(fname=p.input.pdb_in, selection=sel)

###########################################################################
#                              norm occ                                   #
###########################################################################

    if p.params.occ_model == True:
        l.process_message('Normalizing occupancies based on number of models...')
        model_num = len(pdb_object.hierarchy.models())
        new_occ = 1.0 / model_num
        l.show_info('Setting occupancy to {0:.2}'.format(new_occ))
        for at in pdb_object.hierarchy.atoms():
            at.occ = new_occ
    pdb_object.write_hierarchy_to_pdb(output_hierarchy=pdb_object.hierarchy,
                                      out_name=p.output.pdb_out)

###########################################################################
#                              set occ                                    #
###########################################################################

    if p.params.set_occ == True:
        l.process_message('Setting occupancies to input value...')
        new_occ = p.params.occ_val
        l.show_info('Setting occupancy to {0:.2}'.format(new_occ))
        for at in pdb_object.hierarchy.atoms():
            at.occ = new_occ


###########################################################################
#                              set B                                      #
###########################################################################

    if p.params.set_B == True:
        pdb_object.set_B_zero()
        

    pdb_object.write_hierarchy_to_pdb(output_hierarchy=pdb_object.hierarchy,
                                      out_name=p.output.pdb_out)


    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
