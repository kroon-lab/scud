#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

class Plot(object):

    def __init__(self):
        ### Global parameters:
        self.line_data = []
        ### MPL init params:
        mpl.rcParams['axes.labelsize'] = 13
        mpl.rcParams['font.size'] = 13
        mpl.rcParams['legend.fontsize'] = 13
        mpl.rcParams['legend.loc'] = 'best'

    def collect_line_data(self,xlabel,ylabel,label,color,x,y,fname,lw=3,ls='-'):
        '''
        Constructs a dictionairy for line plot data
        '''
        line = {}
        line["xlabel"] = xlabel
        line["ylabel"] = ylabel
        line["label"] = label
        line["color"] = color
        line["x"] = list(x)
        line["y"] = list(y)
        line["fname"] = fname
        line["lw"] = lw
        line["ls"] = ls
        self.line_data.append(line)
        
    def contour2D(self,slice2D,plotfile_name=None,vmin=None,vmax=None):
        '''
        Plots 2D contour of a np.array
        '''
        if plotfile_name == None: plotfile_name = 'plt.eps'
        # Start figure and Axis
        fig = plt.figure(figsize=(13,10))
        ax1 = plt.Axes(fig,[0.0,0.0,1.0,1.0])
        # Turn axes (x,y) off
        ax1.set_axis_off()
        fig.add_axes(ax1)
        # Generate leves based on min and max value
        levels = np.linspace(vmin,vmax,num=15,endpoint=True)
        # Create contour plot, contoured add levels defined above
        img = ax1.contourf(slice2D,levels)
        # Set colormap, could make this user input as well
        img.set_cmap('rainbow')
        # Force square plot
        ax1.set_aspect('equal')
        # Add colorbar
        cbar = plt.colorbar(img,shrink=0.85)
        # Save figure as .eps
        fig.savefig(plotfile_name,format='eps',transparent=True)

    def lines_plot(self,plotfile_name=None):
        '''
        line plot method, plots self.line_data
        '''
        if plotfile_name == None: plotfile_name = 'plt.eps'
        fig = plt.figure(figsize=(9,6))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        for line in self.line_data:
            plt.plot(line["x"],line["y"],label=line["label"],lw=line["lw"],ls=line["ls"],color=line["color"])
            self.set_line_vars(ax1)
        ax1.set_xlabel(self.line_data[0]["xlabel"])
        ax1.set_ylabel(self.line_data[0]["ylabel"])
        fig.savefig(plotfile_name,format='eps',transparent=True)

    def set_line_vars(self,ax1):
        '''
        Extra customization for line plots
        '''
        plt.legend(frameon=False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

    def write_line_data(self):
        '''
        Write line data (dictionairy) as .json
        '''
        import json
        for dat in self.line_data:
            json.dumps(dat,sort_keys=True)
            with open(dat["fname"],'w') as outfile:
                json.dump(dat,outfile)
