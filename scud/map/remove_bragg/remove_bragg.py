#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import uctbx
import iotbx.ccp4_map
from cctbx import sgtbx

import phil as remove_bragg_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass

def run(args=None):
    '''
    Remove Bragg peaks form a .map file. Indices of the Bragg are first
    derived based on the super_size. A leak factor can be added to 
    Filter the neighbouring voxels as wel. Filler parameter determines
    what value will be placed in the selected voxels. 
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='remove_bragg_log.txt')
    l.title("map.remove_bragg module")

###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    # use phil to process input
    working_params = remove_bragg_phil.phil_parse(args=args,log=l)
    p = working_params.remove_bragg

###########################################################################
#                              Read Map                                   #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_in))
    map_in = MAPClass(p.input.map_in)

###########################################################################
#                            Print Info                                   #
###########################################################################

    l.warning('map in info: ')
    map_in.map_info(l=l)

###########################################################################
#                            generating indices                           #
###########################################################################

    l.process_message('Generating Bragg voxel indices...')
    # Map properties needed for bragg index generation
    grid = map_in.maps.unit_cell_grid
    uc = map_in.uc
    offset = ((np.array(grid)-1)/2) + 1
    # Processing input params:
    super_size = p.params.super_size
    leak_num = p.params.leak_num
    if not p.params.median:
        l.process_message('Setting all Bragg voxels to filler value...')
        # heavy for-loop-inception to generate al indices to be masked
        coordlist = make_bragg_list(grid,super_size,leak_num,offset)

###########################################################################
#                      Looping over data array                            #
###########################################################################

        l.process_message('Looping over data array...')
        filler = p.params.filler
        for crd in coordlist:
            map_in.maps.data[crd] = filler

###########################################################################
#                    Applying median filter to map                        #
###########################################################################

    else:
        l.process_message('Applying median filter...')
        median_filter(grid,
                      super_size,
                      leak_num,
                      offset,
                      map_in,
                      l)

###########################################################################
#                      writing new filtered map                           #
###########################################################################

    l.process_message('Writing new map...')
    # Converting data array to double array for output
    dat = map_in.maps.data.as_double()
    iotbx.ccp4_map.write_ccp4_map(
        file_name=p.output.map_out,
        unit_cell = uctbx.unit_cell(parameters = uc),
        space_group=sgtbx.space_group_info("P1").group(),
        map_data=dat,
        labels=flex.std_string(["iotbx.ccp4_map.tst"]))

###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()

def median_filter(grid,
                  super_size,
                  leak_num,
                  offset,
                  map_in,
                  l):
    '''
    Finds Bragg voxels, depending on super_size and offset. Leak num determines the amount of voxels per Bragg. Then check around the Bragg voxels to find a median value to place.
    '''
    l.warning('Will be heavy!')
    # init values
    leak_max = leak_num + 1 # better do the addition only 1 time
    offset = offset.astype(int) # one conversion
    median_box = 1
    median_max = median_box + 1
    
###########################################################################
#                        Generate Bragg indices                           #
###########################################################################
    l.process_message('Generating Bragg List...')
    # xyz are the bragg positions
    bragg_list = [np.array([x,y,z]) + offset
                  for x in range(0 - offset[0], grid[0] - offset[0], super_size)
                  for y in range(0 - offset[1], grid[1] - offset[1], super_size)
                  for z in range(0 - offset[2], grid[2] - offset[2], super_size)]
    l.show_info('{} Bragg positions found in map'.format(len(bragg_list)))
    l.process_message('Applying median filter to Bragg boxes...')

###########################################################################
#                        Generate Bragg boxes                             #
###########################################################################

    for xyz in bragg_list:
        # Find bragg box
        x,y,z = xyz[0],xyz[1],xyz[2]
        bragg_box = [[i,j,k]
                     for i in range(x - leak_num, x + leak_max)
                     for j in range(y - leak_num, y + leak_max)
                     for k in range(z - leak_num, z + leak_max)]
        filtered  = []

###########################################################################
#                  Median filter Bragg box voxels                         #
###########################################################################        

        for ijk in bragg_box:
            # check if in grid
            if indices_in_grid(ijk,grid):
                # xyz (bragg) will get the median of the surrounding box
                if set(ijk) != set(xyz):
                    # find ijk box
                    i,j,k = ijk[0],ijk[1],ijk[2]
                    # create box
                    ijk_box = [[m,n,o]
                               for m in range(i - median_box, i + median_max)
                               for n in range(j - median_box, j + median_max)
                               for o in range(k - median_box, k + median_max)]
                    median_vals = []
                    for mno in ijk_box:
                        # check if in grid and not in bragg box
                        if not element_in_list(mno,bragg_box):
                            if indices_in_grid(mno,grid):
                                median_vals.append(map_in.maps.data[mno])
                    # Calculate median for mno box
                    if len(median_vals) > 0:
                        f = np.median(median_vals)
                        # Filtered contains all medians of bragg box except bragg itself
                        filtered.append(f)
                        map_in.maps.data[ijk] = f
                    else:
                        map_in.maps.data[ijk] = 0
        # calculate median for Bragg reflection
        if len(filtered) > 0:
            map_in.maps.data[xyz] = np.median(filtered)
        else:
            map_in.maps.data[xyz] = 0

def make_bragg_list(grid,super_size,leak_num,offset):
    '''
    Return a list of indices that describe the Bragg positions
    and the positions around it.
    '''
    coordlist = [] # will be returned
    leak_max = leak_num + 1 # better do the addition only 1 time
    offset = offset.astype(int) # one conversion
    # xyz are the bragg positions
    for x in range(0 - offset[0], grid[0] - offset[0], super_size):
        for y in range(0 - offset[1], grid[1] - offset[1], super_size):
            for z in range(0 - offset[2], grid[2] - offset[2], super_size):
                # ijk are the neighbours, accounting for leaking intensity
                for i in range(x - leak_num, x + leak_max):
                    for j in range(y - leak_num, y + leak_max):
                        for k in range(z - leak_num, z + leak_max):
                            # convert from hkl to grid indices
                            ijk = np.array([i,j,k])+offset 
                            # Check if the indices are in the grid
                            if (
                                    0 <= ijk[0] < grid[0] and
                                    0 <= ijk[1] < grid[1] and
                                    0 <= ijk[2] < grid[2]
                            ):
                                # Append to coordlist
                                coordlist.append(ijk)
    return coordlist

def indices_in_grid(i,g):
    '''
    Check if a 3 long index is within a grid
    returns True if it is, False if not
    '''
    if (
            0 <= i[0] < g[0] and
            0 <= i[1] < g[1] and
            0 <= i[2] < g[2]
    ):
        return True
    return False
            
def element_in_list(s,l):
    '''
    Check if an element is in a list
    element can be a 3 long vector
    returns True if it is, False if not
    '''
    s = set(s)
    for a in l:
        if s == set(a):
            return True
    return False
