from __future__ import division
import iotbx.pdb

import numpy as np
from random import randint
import os,sys

import phil as ens2b_phil
from scud.pdb.PDB import PDBClass
from scud.general.log_writer import Log
from scud.general.Plot import Plot

def run(args=None,l=None):
    '''
    Read ensemble (PDB format) calculate average structure and output B-Factor
    calculated from the RMSF
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################

    # init log file
    if l == None:
        l = Log(log_name='ens2b_log.txt')
    l.title('ens2b module')

###########################################################################
#                      Process input Params                               #
###########################################################################

# use phil to process input
    working_params = ens2b_phil.phil_parse(args=args,log=l)
    p = working_params.ens2b

###########################################################################
#                      Read and analyse pdb_in                            #
###########################################################################

    l.process_message('Reading and analyzing PDB file / Structure...')
    pdb_object = PDBClass(fname=p.input.pdb_in, selection="not (resname HOH) and not (resname CL) and not (resname NA)")

###########################################################################
#                        Calc RMSF from ensemble                          #
###########################################################################

    x,y = pdb_object.ensemble_to_B_factor(ensemble_hierarchy=pdb_object.hierarchy)

###########################################################################
#                         Plot RMSF as B-fact                             #
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
    l.show_info('plot file name: {}'.format(p.output.plot_out))
    pt.lines_plot(plotfile_name=p.output.plot_out)
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    return l
