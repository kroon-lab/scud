#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os
import json

from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import uctbx
import iotbx.ccp4_map
from cctbx import sgtbx

import phil as radial_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass

def run(args=None):
    '''
    Calculate radial averages of maps, and subtract them
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='radial_log.txt')
    l.title("map.radial module")

###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    # use phil to process input
    working_params = radial_phil.phil_parse(args=args,log=l)
    p = working_params.radial

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
#          Find all d2's and intensities needed for radial profile        #
###########################################################################

    if p.input.radial_in == None:
        l.process_message('no radial profile given...')
        # Generate a list with only non bragg voxels (+- leak number)
        d2_i_list = make_d2_i_list(grd = map_in.grd,
                                   super_size = p.params.super_size,
                                   leak_num = p.params.leak_num,
                                   l = l,
                                   map_in = map_in)
        l.show_info('d2 list of {} elements made'.format(len(d2_i_list)))

###########################################################################
#                      Construct radial profile                           #
###########################################################################

        # generate radial average
        num_bins = 50 # Make input param (phil)!!!!!!
        radial = make_radial(l = l,
                             num_bins = num_bins,
                             d2_i_list = d2_i_list)
    else:
        with open(p.input.radial_in,'r') as radial_file:
            radial = json.load(radial_file)

###########################################################################
#                          Plot radial profile                            #
###########################################################################
            
    plot_radial(l = l,
                radial = radial)

###########################################################################
#               Subtract radial profile and save new map                  #
###########################################################################

    map_in.subtract_radial(l = l,
                           radial = radial) 
    
##########################y#################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()

def plot_radial(l = None,
                radial = None,
                outfile = 'radial.eps'):
    '''
    Plot radial profile using scud.plt.Plot
    '''
    from scud.plt.Plot import Plot
    l.process_message("Plotting radial profile")
    pt = Plot()
    pt.collect_line_data("resolution ($A^{-2}$)",
                         "Intensity",
                         "radial average (map)",
                         "blue",
                         radial["d2_mean"],
                         radial["i_mean"],
                         outfile)
    pt.lines_plot()

def make_radial(l = None,
                num_bins = 20,
                d2_i_list = None):
    '''
    Gets a list with reso's (d^2) and intensities, bins list and returns
    a radial average list (binned list) containing bin limits
    '''
    l.process_message('generating radial average')
    # Sort list
    d2_i_list = sorted(d2_i_list, key=lambda x: x[0])
    l.show_info('d2_i_list sorted\n')

    # setup binning:
    binned = [] # Will be returned containing stats
    d2_range = d2_i_list[-1][0] - d2_i_list[0][0] # Needed for binning
    d2_bin_size = d2_range / num_bins # d2 range per bin
    d2_min = d2_i_list[0][0] # first minimal value

    # Actual binning
    for n in range(num_bins):
        # find bin values:
        bin_l = [i for i in d2_i_list if d2_min < i[0] < d2_min + d2_bin_size ] 
        d2_min = bin_l[-1][0] # set next bin start
        d2,i = zip(*bin_l) # separate d2 and i
        mean_i = np.mean(i) # Average bin intensity
        mean_d2 = np.mean(d2) # average resolution (for plotting)
        d2_range = [np.min(d2),np.max(d2)] # bin limits
        binned.append([mean_d2,d2_range,mean_i])
        if n % 5 == 0: l.show_info('bin {} finished'.format(n+1))
        
    # Save radial profile info in .json format
    radial = json_dump_radial(binned = binned,
                              l = l)
    
    # Return results
    return radial

def json_dump_radial(binned = None,
                     l = None):
    """
    Dump radial profile as json file
    """
    l.process_message('writing radial profile to json file')
    # Final dictionary
    radial = {}

    # Fill dict
    zp_d2_mean, zp_d2_range, zp_i_mean = zip(*binned)
    radial["d2_mean"] = zp_d2_mean
    radial["d2_range"] = zp_d2_range
    radial["i_mean"] = zp_i_mean

    # Make json format
    json.dumps(radial,sort_keys=True)

    # save to json
    with open("radial.json",'w') as outfile:
        json.dump(radial,outfile)
    l.show_info('Radial profile info written to radial.json')

    # Return dict
    return radial
    
def make_d2_i_list(grd = None,
                   super_size = None,
                   leak_num = None,
                   l = None,
                   map_in = None):
    '''
    Generate an index list of non-Bragg / non-leaked indices of the map
    check if these indices are within max reso, and save reso and intensity
    to list
    '''
    l.process_message('Generate index list...')

    # Shift from map 0,0,0 to reciprocal space 0,0,0
    origin = map_in.origin_shift

    # center coordinate in reciprocal space (d^2)
    center = np.array([0,0,0])

    # data from map
    dat = map_in.maps.data
    
    # max resolution
    max_d2 = map_in.max_reso_d2()
    l.show_info('Max reso (a^-2): {}'.format(max_d2))
    
    # Account for different lengths of map, and origin shift
    limit_list = range(-leak_num,leak_num+1)
    limit_list_x = [(x + origin[0]) % super_size for x in limit_list]
    limit_list_y = [(y + origin[1]) % super_size for y in limit_list]
    limit_list_z = [(z + origin[2]) % super_size for z in limit_list]
    
    # Results are written into
    d2_i_list = []
    
    # Huge loop, generate all indices possible, check if they are Bragg
    # check if the resolution is not exceeding max resolution 
    # Then save reso and intensity to final list
    for x in range(0, grd[0]):
        if x % 50 == 0: l.show_info('x = {}'.format(x))
        if x % super_size not in limit_list_x:
            tmp = list(flex.nested_loop((x,0,0),(x,grd[1],grd[2])))
            for ind in tmp:
                ind = np.array(ind)
                d2 = map_in.ind_to_d_2(center,ind)
                # append if under max reso
                if d2 < max_d2: d2_i_list.append([d2,dat[ind]])
            continue
        for y in range(0, grd[1]):
            if y % super_size not in limit_list_y:
                tmp = list(flex.nested_loop((x,y,0),(x,y,grd[2])))
                for ind in tmp:
                    ind = np.array(ind)
                    d2 = map_in.ind_to_d_2(center,ind)
                    # append if under max reso
                    if d2 < max_d2: d2_i_list.append([d2,dat[ind]])
                continue            
            for z in range(0, grd[2]):
                if z % super_size not in limit_list_z:
                    # calculate resolution
                    ind = np.array([x,y,z])
                    d2 = map_in.ind_to_d_2(center,ind)
                    # append if under max reso
                    if d2 < max_d2: d2_i_list.append([d2,dat[ind]])
    return np.array(d2_i_list)
    
def do_statistics(ls): # Needs new name - OLD
    '''
    Do statistics on map values based on resolution
    later should be used for generating spherical radial average
    '''
    # Sort list
    ls = sorted(ls, key=lambda x: x[0])
    ls_binned = [] # Will be returned containing stats
    n_bins = 50 # Should be made user input
    d2_range = ls[-1][0] - ls[0][0] # Needed for binning
    d2_bin_size = d2_range / n_bins # d2 range per bin
    d2_min = ls[0][0] # first minimal value
    for n in range(n_bins):
        bin_l = [i for i in ls if d2_min < i[0] < d2_min + d2_bin_size ] # find bin values
        d2_min = bin_l[-1][0] # set next bin start
        d2,i = zip(*bin_l) # separate d2 and i
        # Rest not really needed
        mean_d2 = np.mean(d2)
        reso = 1./np.sqrt(mean_d2)
        mean_i = np.mean(i)
        std = np.std(i)
        d2_range = [np.min(d2),np.max(d2)]
        min_reso = 1./np.sqrt(ls[0][0])
        ls_binned.append(['%.2f'%mean_d2, '%.2f'%reso, mean_i, std,i,'%.2f'%min_reso,d2_range]) # Return labels as string
    return ls_binned # Return all stats      
    
