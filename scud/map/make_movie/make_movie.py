#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys,os

import phil as make_movie_phil
from scud.general.log_writer import Log
from scud.map.MAP import MAPClass
from scud.map.MAP import MapMapClass
from scud.plt.Plot import Plot
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib

def run(args=None):
    '''
    Make a movie odf slicing through 2 maps
    '''
###########################################################################
#                         Start log file                                  #
###########################################################################
    # init log file
    l = Log(log_name='make_movie_log.txt')
    l.title("map.make_movie module")
    
###########################################################################
#                      Process input Params                               #
###########################################################################

    l.process_message('Processing input...\n')
    working_params = make_movie_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.make_movie

###########################################################################
#                              Read Map                                   #
###########################################################################

    l.process_message('Opening map {}'.format(p.input.map_1))
    map_1 = MAPClass(p.input.map_1)
    l.process_message('Opening map {}'.format(p.input.map_2))
    map_2 = MAPClass(p.input.map_2)

###########################################################################
#                            Print Info                                   #
###########################################################################
    # Show map info
    map_1.map_info(l=l)
    map_2.map_info(l=l)

###########################################################################
#                            PREP Data                                    #
###########################################################################

    # Load MapMapClass
    mm = MapMapClass(map_1 = map_1,
                     map_1_dat = map_1.maps.data,
                     map_2 = map_2,
                     map_2_dat = map_2.maps.data,
                     p = p,
                     l = l)
    # Resize Maps
    mm.resize_maps(p = p,
                   l = l)

    msk = (mm.dat_1 == -1000)
    dat_1 = mm.dat_1.set_selected(msk, 0)
    dat_2 = mm.dat_2.set_selected(msk, 0)    
    
    # convert to numpy array
    dat_1 = dat_1.as_numpy_array()
    dat_2 = dat_2.as_numpy_array()
    
###########################################################################
#                        Prep tmp folder                                  #
###########################################################################

    # Make tmp folder for separate images
    tmp_folder = 'tmp_movie'
    if os.path.isdir(tmp_folder):
        os.rmdir(tmp_folder)
    os.mkdir(tmp_folder)

###########################################################################
#                            Plot Slices                                  #
###########################################################################

    # Save file names
    flist = []
    # Loop over l-slices
    for i in range(np.shape(dat_1)[2]):
        # Sequencial file names
        fname = tmp_folder+'/{:0>3}.png'.format(i)
        flist.append(fname)
        # If log or not take log, rot 90 deg for proper orientatino
        if p.params.log:
            t_1 = np.log(np.rot90(dat_1[:,:,i]))
            t_2 = np.log(np.rot90(dat_2[:,:,i]))
        else:
            t_1 = np.rot90(dat_1[:,:,i])
            t_2 = np.rot90(dat_2[:,:,i])

        # Prep plot
        fig = plt.figure()
        # 2 subplots
        gs = gridspec.GridSpec(1,2)

        # Plot 1
        ax1 = fig.add_subplot(gs[0,0])
        ax1.set_axis_off()
        ax1.imshow(t_1,cmap='Greys')
        ax1.set_aspect('equal')

        # Plot 2
        ax2 = fig.add_subplot(gs[0,1])
        ax2.imshow(t_2,cmap='Greys')
        ax2.set_axis_off()
        ax2.set_aspect('equal')

        # Layour info
        plt.tight_layout()
        # Save fig
        fig.savefig(fname,format='png',transparent=True)
    
###########################################################################
#                            Create GIF                                   #
###########################################################################

    movie_name = p.output.movie_out # Should be mp4!!!!
    # Make movie
    os.system('ffmpeg -i {}/%03d.png {}'.format(tmp_folder,movie_name))
    
###########################################################################
#                            Close Log File                               #
###########################################################################

    l.close_log()
