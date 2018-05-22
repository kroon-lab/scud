#!/usr/bin/env cctbx.python

from __future__ import division
import iotbx
from iotbx import pdb
import numpy as np
from scitbx.array_family import flex

from scud.general.method_library import vecDistSquare


class ERClass(object):
    '''
    TODO: 
    * Input: hierarchy or filename
    * Method for sc reduction
    * Method to create test structure factors from ensemble descibing motion and prepare input supercell forensemble refinement
    * Method to create a ER-parameter file needed to run ensemble refinement with outputted structure
    '''
    def __init__(self):


    def reduce_supercell(self):
        '''
        Reduce supercell. Use uctbx.fractionalize to fractionalize coordinates to the 'original' unit cell. Reduce coordinates of center of mass of molecules to be as close as possible to 0-1 range. Then orthogonalize, write to file
        extra feature: reduce to asu?

        input:
        * hierarchy of supercell
        * supercell size (integer or 3-long integer list)
        * extra: space group symmetry to reduce to original ASU
        '''

    def data_from_enemble(self):
