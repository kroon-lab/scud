#!/usr/bin/env cctbx.python

from __future__ import division
from iotbx import mtz
#import iotbx.ccp4_map
from cctbx import miller
#from scitbx.array_family import flex
from cctbx.array_family import flex
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import crystal
import numpy as np
#import makeMap
import math
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.mlab import griddata
#from libtbx import easy_pickle
import phil as ewald_phil
from scud.general.log_writer import Log
from scud.mtz.MTZ import MTZClass
#import os.path

def Ewaldcenter(center):
    if center != None:
#    if os.path.isfile("beam.in") == True:
#        f = open("beam.in")
#        Loes = f.readline().split()
        center = center.split(',')
        center = map(float,center)
#        f.close()
#        Loes = float(Loes)
    else:
        center = (1.0,0.0,0.0)
    Loes = (0.542114, -0.261678, -0.839945)
    Tim = (0.148171,1.00058,0.211549)
    Perez = (-0.319537,0.21346,1.04254)
    return center

def wavelength():
    ESRF = 0.9677
    Perez = 0.9
    return ESRF

def distance_3D(vec1,vec2):
    '''
    return 3D distance vector from 2 vectors (xyz positions)
'''
    return math.sqrt( (vec1[0] - vec2[0])**2 + (vec1[1] - vec2[1])**2 + (vec1[2] - vec2[2])**2 )



def run(args=None, l=None):

    '''
    Convert structure factor file (mtz) to ccpp4 type map.
    '''

#    help(iotbx.mtz.object())
#    quit()
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='ewald_log.txt')
    l.title("ewald module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = ewald_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.ewald
    
    '''
    Extract intensities that lie on the Ewald sphere and project them on a plane (detecor)
    '''
    '''
    For array in mtz check if Miller index is on ewald sphere
    '''
    # Read MTZ file and extract info + miller index

    mtz_in = MTZClass(p.input.mtz_in,p.params.array)
    print 'Space group',mtz_in.symmetry.space_group_info()
    print 'Unit cell',mtz_in.uc
    ma = mtz_in.miller_array
    diff_miller = ma
    # Determines ewald sphere:
    center  = p.params.center
    center = np.array(Ewaldcenter(center))
    print 'center',center
    limit = 0.001
    limit = 0.003
    refDist= 1./wavelength()
    hkl_I = []
    # Loop over all SFs
    print 'Looping over all SFs'
    for i,sf in enumerate(diff_miller): 
        # Factional coordinates from Miller index
        if i%1000000 == 0:
               li = len(str(i))
               print i
        frac_hkl = mtz_in.uc.reciprocal_space_vector(sf[0])
        # Distance from ewald sphere center
        distance = distance_3D(frac_hkl,center)
        # Check if within shell of Ewald sphere
        ind = np.array([sf[0][0],sf[0][1],sf[0][2]])
        if distance < (refDist + limit) and distance > (refDist - limit):
            hkl_I.append([tuple(ind),sf[1]+100])
        else:
            hkl_I.append([tuple(ind),0.0])
        if sf[0][2] != 0 :
            centro_hkl = mtz_in.uc.reciprocal_space_vector((-1*sf[0][0],-1*sf[0][1],-1*sf[0][2]))            # Generate centro symmetrical miller
            ind = np.array([-1*sf[0][0],-1*sf[0][1],-1*sf[0][2]])
            distance = distance_3D(centro_hkl,center)
            if distance < (refDist + limit) and distance > (refDist - limit):
                hkl_I.append([tuple(ind),sf[1] + 100.0])
            else:
                hkl_I.append([tuple(ind),0.0])

    mi = flex.miller_index(zip(*hkl_I)[0])
    ms = miller.set(
        crystal_symmetry = mtz_in.symmetry,
        anomalous_flag = True,
        indices = mi)
    new_data = flex.double(zip(*hkl_I)[1])
    ma = miller.array(ms,
                      data=new_data)
    mtz_dataset = ma.as_mtz_dataset(column_root_label='IDFF',column_types='J')
    mtz_dataset.mtz_object().write('image.mtz')
    
                





