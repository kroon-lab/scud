#!/usr/bin/env cctbx.python

from __future__ import division
import os, sys
import numpy as np

import iotbx
from iotbx import pdb
import numpy as np
from scitbx.array_family import flex

from scud.pdb.supercell import supercell
from scud.mtz.mtz2map import mtz2map

class EXPClass(object):
    '''
    TODO: 
    * docstrings :P
    * ensemble trimmer based on either:
        * Structure Factor variance (example: lowest 10% R-value structures with respect to the average structure?)
        * RMSD, X structures with lowest RMSD compared to average structure?
    '''
    
    def __init__(self):
        t = 10

    def calc_diffuse_and_map(self,
                             pdb = None,
                             supercell_num = None,
                             size_h = None,
                             size_k = None,
                             size_l = None,
                             Ncpu = None,
                             write_pdb = None,
                             l = None):
        '''
        Calculate diffuse scattering using scud.pdb.supercell
        rename output a nd remove temporary files
        Convert mtz IDFF column to map
        '''
        # Input for diffuse scatterign calculation
        supercell_arg = ['pdb_in={}'.format(pdb),
                         'supercell_num={}'.format(supercell_num),
                         'size_h={}'.format(size_h),
                         'size_k={}'.format(size_k),
                         'size_l={}'.format(size_l),
                         'Ncpu={}'.format(Ncpu),
                         'write_pdb={}'.format(write_pdb)]
        supercell.run(args=supercell_arg,
                      l = l)
        # Remove tmp files rename .mtz
        os.system('rm tmp_*.mtz')
        os.system('mv result.mtz {}.mtz'.format(pdb[:-4]))
        # Input for mtz2map calculation
        mtz2map_arg = ['mtz_in={}.mtz'.format(pdb[:-4]),
                       'array=IDFF',
                       'map_out={}_IDFF.map'.format(pdb[:-4])]
        mtz2map.run(args = mtz2map_arg,
                    l = l)
