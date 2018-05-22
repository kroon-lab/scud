from __future__ import division
import iotbx.pdb

import numpy as np
from random import randint
import os,sys

import phil as b_phil
from scud.pdb.PDB import PDBClass
from scud.general.log_writer import Log
from scud.general.Plot import Plot

def run(args):
    '''
    Read ensemble (PDB format) calculate average structure and output B-Factor
    calculated from the RMSF
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    l = Log(log_name='b_log.txt')
    l.title('b module')

###########################################################################
#                      Process input Params                               #
###########################################################################

# use phil to process input
    working_params = b_phil.phil_parse(args=args,log=l)
    p = working_params.b

###########################################################################
#                      Read and analyse pdb_in                            #
###########################################################################

    l.process_message('Reading and analyzing PDB file / Structure...')
    if p.params.av_all_res == False:
        sel = '(name CA) not (resname HOH) and not (resname CL) and not (resname NA) and (altid " " or altid "A")'
    else:
        sel = 'not (resname HOH) and not (resname CL) and not (resname NA) and (altid " " or altid "A")'
    pdb_object = PDBClass(fname=p.input.pdb_in, selection=sel)

###########################################################################
#                        Extract B facts                                  #
###########################################################################

    # Get residue numbers
    x = [int(i.resseq) for i in pdb_object.hierarchy.models()[0].chains()[0].residues()]
    if p.params.av_all_res==False:
        y = [float(i.b) for i in pdb_object.hierarchy.models()[0].chains()[0].atoms()]
        ylab = 'Refined B (CA)'
        fn = 'b_ca.json'
    else:
        y = []
        for res in pdb_object.hierarchy.models()[0].chains()[0].residues():
            y.append(np.mean([a.b for a in res.atoms()]))
        ylab = 'Refined B (averaged)'
        fn = 'b_av.json'
            
###########################################################################
#                         Plot B-fact                                     #
###########################################################################

    l.process_message('Plotting B-Factors...')
    #### Plot
    pt = Plot()
    pt.collect_line_data(xlabel='Residue Number',
                         ylabel='B-factor ($\AA^{2}$)',
                         label=p.output.plot_dat[:-4],
                         color=(0.8,0.4,0.0),
                         x=x,
                         y=y,
                         fname=p.output.plot_dat)
    pt.write_line_data()
    l.show_info('plot file name: {}'.format(p.output.plot_out[:-4]+fn))
    pt.lines_plot(plotfile_name=p.output.plot_out)
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
