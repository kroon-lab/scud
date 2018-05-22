#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import mmtbx.utils
from libtbx import easy_mp
from iotbx import mtz
from cctbx import maptbx

from scitbx.array_family import flex

import phil as map_radial_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass
from scud.map.MAP import MapMapClass
from scud.plt.Plot import Plot

def run(args=None, l= None):
    '''
    Check additivity of motion in reciprocal space, subtract mtz files
    '''
    
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    if l == None:
        l = Log(log_name='map_radial_log.txt')
    l.title("exp.map_radial module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = map_radial_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.map_radial

    map_in = MAPClass(p.input.map_in)
    map_in.map_info(l = l)
    map_dat = map_in.maps.data

###########################################################################
#                           Radial Operarations                           #
###########################################################################

    l.process_message('Constructing radial profiles...')
    # Prep bin sizes based on resolution and number of bins
    d = map_in.reso[0] # Max resolution
    l.show_info('Max resolution: {:.2f}'.format(d))
    n = p.params.n_bins # Number of bins
    # resolution list to bin in (1/d^3):
    d_pwr = np.power(list(np.arange(0,(1/d**3),(1/d**3)/(1.0*n))) + [1/d**3],(1/3.))
    
    x_axis = (d_pwr[1:] + d_pwr[:-1])/2. # X-axis used for plotting (later)
    
    # Different formats of grid needed, 1d, 3d and flex.grid
    grd = map_in.grd
    grd_1d = np.product(grd)
    grd_reshape = flex.grid(grd)
    
    # Needed for gridding around sites:
    unit_cell = map_in.unit_cell
    uc = map_in.uc
    # Not measured values
    ignore_val = 0
    invalid_mask = (map_dat == ignore_val)
    
    # Initialise blank mask (False)
    mask_lower_binary = flex.bool(grd_1d, False)
    mask_lower_binary.reshape(grd_reshape)
    
    # Will be updated with new lower level
    d_star_min = 0.0
    
    # Containers
    m1_I_l, m2_I_l, I_lists = [],[],[]
    reso_CC_ls = []
    stat_string_ls = []
    
    # For subtraction and output:
    radial_map_1 = map_dat.deep_copy()
    
    # Will become aniso maps:
    m1 = map_dat
    
    # Loop over resolution limits
    for i,d_reso in enumerate(d_pwr[1:]):
        # Upper limits
        mask_upper = maptbx.grid_indices_around_sites(unit_cell  = unit_cell,                                           # uctbx unit_cell object
                                                      fft_n_real = grd,                                                 # tuple (array size)
                                                      fft_m_real = grd,                                                 # tuple
                                                      sites_cart = flex.vec3_double([np.array(uc)[0:3]/2.]),            # origin of array (flex.vec3.double)
                                                      site_radii = flex.double(1,d_reso))                               # Reso of shell (radius around sites_cart)
        
        # Convert mask to binary
        mask_upper_binary = flex.bool(grd_1d, False).set_selected(mask_upper, True)
        mask_upper_binary.reshape(grd_reshape)
        # Masking everythin not between lower and upper limit
        radial_mask = mask_upper_binary.deep_copy().set_selected(mask_lower_binary,False)
        # Add invalid mask
        total_mask = radial_mask.deep_copy().set_selected(invalid_mask,False)
        
        # Select all values in resolution shell 
        m1_bin_I = m1.select(total_mask.as_1d())
        # Calculate mean I
        m1_I_mean = flex.mean(m1_bin_I)
        
        stat_s = '{:.2f}-{:.2f} ({}): {:.2f}'.format(d_star_min,d_reso,m1_bin_I.size(),m1_I_mean)
        l.show_info(stat_s)
        stat_string_ls.append(stat_s)
        m1_I_l.append(m1_I_mean)
        
        # Subtract radial average from masked part of aniso_map:
        radial_map_1.set_selected(total_mask,m1_I_mean)
        
        # Move upper mask to lower mask of next shell
        mask_lower_binary = mask_upper_binary
        # Update lower resolution limit
        d_star_min = d_reso
        
    # Append to I_list container
    print m1_I_l

    # Subtract radial mask from dat set invalid to 0 and write map
    l.process_message('Creating and writing aniso map...')
        
    # Create aniso maps
    #if not p.params.map_1_exp:
    m1 = m1 - radial_map_1
    
    # Final aniso maps invalid mask set to 0 because intensities are distributed around 0:
    m1.set_selected(invalid_mask,0)

    map_in.write_map(dat = m1,
                     name = p.input.map_in[:-4]+'_aniso.map',
                     l = l)

    # Write radial trace to file
    pt = Plot()
    pt.collect_line_data(xlabel='Reso (1/d$^{3}$)',
                         ylabel='Intensity',
                         label=p.input.map_in[:-4],
                         color=(0.2,0.2,0.2),
                         x=x_axis,
                         y=m1_I_l,
                         fname=p.input.map_in[:-4]+'_radial_line.json')
    pt.write_line_data()

    with open(p.input.map_in[:-4]+'_bins.txt','w') as f:
        for l in stat_string_ls:
            print >> f, l
    
    return l

