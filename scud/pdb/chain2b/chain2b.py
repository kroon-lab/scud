from __future__ import division
import iotbx.pdb
import iotbx.cif

import numpy as np
from random import randint
import os,sys

import phil as chain2b_phil
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
    l = Log(log_name='chain2b_log.txt')
    l.title('chain2b module')

###########################################################################
#                      Process input Params                               #
###########################################################################

# use phil to process input
    working_params = chain2b_phil.phil_parse(args=args,log=l)
    p = working_params.chain2b

###########################################################################
#                          Read  pdb_in                                   #
###########################################################################

    if p.input.cif_in != None:
        struct = iotbx.pdb.input(p.input.cif_in)
        l.process_message('CIF file read...')
    if p.input.pdb_in != None:
        struct = iotbx.pdb.input(file_name=p.input.pdb_in)
        l.process_message('PDB file read...')
    l.process_message('Converting to hierarchy...')
    symmetry= struct.crystal_symmetry()
    hierarchy = struct.construct_hierarchy()

    # Replacing B-factors by 1/i i=chain number
    if not p.params.per_model:
        for i,chain in enumerate(hierarchy.chains()):
            new_b = i / 10.
            for at in chain.atoms():
                at.b = new_b

    # Replace B-factor based on model number
    if p.params.per_model:
        cl = iotbx.pdb.systematic_chain_ids()
        l.process_message('Replacing B-factor column per model...')
        for i, mdl in enumerate(hierarchy.models()):
            new_b = float(i)
            for ch in mdl.chains():
                ch.id = cl[i]
            for at in mdl.atoms():
                at.b = new_b

    l.process_message('B-factors altered...')
    # Writing to file
    if p.output.pdb_out != None:
        l.process_message('Writing to PDB file...')
        hierarchy.write_pdb_file(file_name=p.output.pdb_out,
                                 crystal_symmetry=symmetry,
                                 anisou=False)     
    if p.output.cif_out != None:
        l.process_message('Writing to CIF file...')
        from iotbx.cif import model
        cif = model.cif()
        cif_input = hierarchy.as_cif_input(crystal_symmetry=symmetry)
        cif_block = cif_input.cif_block
        cif["supercell"] = cif_block
        cif.show(out=open(p.output.cif_out,'w'))
    
###########################################################################
#                            Close Log File                               #
###########################################################################
    l.process_message('Done...')
    l.close_log()
