#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

class Plot(object):

    def __init__(self,dark=False):
        ### Global parameters:
        self.line_data = []
        self.bar_data = []
        ### MPL init params:
        if dark:
            plt.style.use('dark_background')
        else:
            plt.style.use('bmh')
        mpl.rcParams['axes.labelsize'] = 13
        mpl.rcParams['font.size'] = 13
        mpl.rcParams['legend.fontsize'] = 13
        mpl.rcParams['legend.loc'] = 'best'

    def bar_plot(self,
                 plotfile_name=None,
                 show=False):
        '''
        Plot bar
        '''
        if plotfile_name == None: plotfile_name = 'plt.eps'
        fig = plt.figure(figsize=(9,6))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        for bar in self.bar_data:
            plt.bar(bar["x"],bar["y"],label=bar["label"],width=bar["width"],color=bar["color"])
#            self.set_line_vars(ax1)
        ax1.set_xlabel(self.bar_data[0]["xlabel"])
        ax1.set_ylabel(self.bar_data[0]["ylabel"])
        fig.savefig(plotfile_name,format='eps',transparent=False)
        if show: plt.show()

    def collect_bar_data(self,xlabel,ylabel,label,color,x,y,width,fname):
        '''
        Constructs a dictionairy for bar plot data
        '''
        bar = {}
        bar["xlabel"] = xlabel
        bar["ylabel"] = ylabel
        bar["label"] = label
        bar["color"] = color
        bar["x"] = list(x)
        bar["y"] = list(y)
        bar["width"] = list(width)
        bar["fname"] = fname
        self.bar_data.append(bar)
        
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

    def colorbar(self,
                 nm = None,
                 vmin = None,
                 vmax = None,
                 lb = 'Intensity'):
        fig = plt.figure(figsize=(1,5))
        #axes [left,bottom,width,height]
        ax = fig.add_axes([0.15,0.15,0.1,0.7])
        cmap = 'Greys'
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        cb = mpl.colorbar.ColorbarBase(ax,
                                       cmap = cmap,
                                       norm = norm,
                                       orientation='vertical')
        cb.set_label(lb)
        fig.savefig(nm,format='eps',transparent=True, bbox_inches='tight')
        
    def contour2D(self,slice2D,
                  plotfile_name=None,
                  vmin=None,
                  vmax=None,
                  cnt_num=None,
                  png_format=False):
        '''
        Plots 2D contour of a np.array
        '''
        if plotfile_name == None: plotfile_name = 'plt.eps'
        if cnt_num == None: cnt_num=15
        show_type = 'imshow'
        # Start figure and Axis
        fig = plt.figure(figsize=(9,9))
        fig.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0,wspace=0.0,hspace=0.0)
        ax1 = fig.add_subplot(111)
        # Turn axes (x,y) off
        ax1.set_axis_off()
        # Start plot
        # Generate leves based on min and max value
        if show_type == 'contour':
            levels = np.linspace(vmin,vmax,num=cnt_num,endpoint=True)
            # Create contour plot, contoured add levels defined above
            img = ax1.contourf(slice2D,levels)
            # Set colormap, could make this user input as well
        if show_type == 'imshow':
            img = ax1.imshow(np.rot90(slice2D), vmin=vmin, vmax=vmax)
        img.set_cmap('Greys')
        # Force square plot
        ax1.set_aspect('equal')
        # Add colorbar
        colorbar = False
        if colorbar:
            cbar = plt.colorbar(img,shrink=0.85)
        if png_format:
            fig.savefig(plotfile_name,format='png',transparent=True)
        else:
            # Save figure as .eps
            fig.savefig(plotfile_name,format='eps',transparent=True, bbox_inches='tight')

    def histogram(self,
                  outname='histogram.eps',
                  dat = None,
                  bin_num= 50):
        '''
        Plot histogram of a 1D list of numbers.
        '''
        fig = plt.figure(figsize=(9,6))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        # Histo plot
        n, bins, patches = plt.hist(dat, bin_num, normed=1)
        ax1.set_xlabel('Value')
        ax1.set_ylabel('Density')
        fig.savefig(outname,format='eps',transparent=False)

    def lines_plot(self,
                   plotfile_name=None,
                   show=False):
        '''
        line plot method, plots self.line_data
        '''
        if plotfile_name == None: plotfile_name = 'plt.eps'
        fig = plt.figure(figsize=(9,6))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        for line in self.line_data:
            plt.plot(line["x"],line["y"],label=line["label"],ls=line["ls"])
            self.set_line_vars(ax1)
        ax1.set_xlabel(self.line_data[0]["xlabel"])
        ax1.set_ylabel(self.line_data[0]["ylabel"])
        fig.savefig(plotfile_name,format='eps',transparent=False)
        if show: plt.show()
        
    def projection(self,
                   dat = None,
                   outname = 'projection.eps',
                   diff = True,
                   l = None):
        '''
        Plot projection of Map on one plane. Find largest value from zero
        use that to color.
        '''
        # Find max value from zero (negative or positive)
        if diff:
            minmax = np.max(np.abs(dat))
        # Create canvas
        fig = plt.figure(figsize=(13,13))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        # No grid
        ax1.grid(b=False)
        dat = np.rot90(dat)
        # Plot projection
        if diff:
            img = ax1.imshow(dat,clim=(-minmax,minmax))
        else:
            img = ax1.imshow(dat,clim=(np.min(dat),np.max(dat)))
        # Color map blue-white-red
        if diff:
            img.set_cmap('bwr')
        else:
            img.set_cmap('bone')
        # Add colorbar
        cbar = plt.colorbar(img,shrink=0.85)
        # THIS ONE DOES NOT WORK
        ax1.set_aspect('equal')
        # No axes ticks and or labels
        ax1.set_axis_off()
        # Save figure
        fig.savefig(outname,format='eps',transparent=True)
        
    def QQ_plot(self,
                dat_1=None,
                dat_2=None,
                xlab='x',
                ylab='y',
                out_name='qq.eps',
                l=None):
        '''
        Calculate slope of 2 data sets, 
        scale dat_1 to dat_2 using slope and intersect
        return slope and intersect
        '''
        # Use percentile (25%)
        perc = 0.25*len(dat_1)
        # polyfit two datasets
        poly = np.polyfit(dat_1[perc:-perc],dat_2[perc:-perc],1)
        # Scale dat_1
        dat_1 = (dat_1*poly[0])+poly[1]
        # Prep figure
        fig = plt.figure(figsize=(9,6))
        fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
        ax1 = fig.add_subplot(111)
        # Scatter plot dat_1 vs dat_@
        plt.scatter(dat_1,
                    dat_2,
                    color='k',
                    marker='.')
        # reference line dat_1 vs dat_1
        plt.plot(dat_1,
                 dat_1,
                 color='r',
                 label='x=y')
        # Set other plot parameters
        self.set_line_vars(ax1)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        # save plot
        fig.savefig(out_name,format='eps',transparent=False)
        # return scaling parameters
        return poly[1],poly[0]
        
    def read_line_data(self,fname=None):
        '''
        open and read json file
        '''
        import json
        with open(fname,'r') as infile:
            d = json.load(infile)
        self.line_data.append(d)

    def set_line_vars(self,ax1):
        '''
        Extra customization for line plots
        '''
        ax1.legend(frameon=False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

    def sub_lines_plot(self,plotfile_name=None):
        '''
        line plot method, plots self.line_data
        '''
        # Cycler to cycle over colors
        from cycler import cycler
        # If no name is given
        if plotfile_name == None: plotfile_name = 'plt.eps'
        # number of sub plots
        n_plots = len(self.line_data)
        # Create figure and subplots
        fig,axarr = plt.subplots(n_plots,sharex=True)
        # Adjust margins
        fig.subplots_adjust(left=0.05,bottom=0.1,right=0.97,top=0.95,wspace=0.1,hspace=0)
        # Create cycles from colors in color theme
        color_cycle = cycler(color=mpl.rcParams['axes.prop_cycle'])
        # counter for recognizing the last plt
        n = 0
        # Iterate over color_cycle and data
        for i,line in zip(color_cycle,self.line_data):
            # Plot line data, use color from theme color params cycle
            axarr[n].plot(line["x"],line["y"],label=str(line["label"]),ls=line["ls"],color=i['color'])
            # Set extra parameters
            self.set_line_vars(axarr[n])
            # For bottom plot, show xticks,the rest is off
            if n != n_plots-1:
                axarr[n].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
            else:
                axarr[n].tick_params(axis='x',which='both',bottom='on',top='off',labelbottom='on')
            # Hide all yticks
            axarr[n].tick_params(axis='y',which='both',left='off',right='off')
            # Count to next plot
            n+=1
        # remove tick labels (numbers along axes)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],visible=False)
        plt.setp([a.get_yticklabels() for a in fig.axes[:]],visible=False)
        # y- and x-labels
        fig.text(0.025,0.5,self.line_data[-1]['ylabel'],ha='center',va='center',rotation='vertical')
        axarr[-1].set_xlabel(self.line_data[-1]['xlabel'])
        # save figure
        fig.savefig(plotfile_name,format='eps',transparent=False)

    def write_bar_data(self):
        '''
        Write line data (dictionairy) as .json
        '''
        import json
        for dat in self.bar_data:
            json.dumps(dat,sort_keys=True)
            with open(dat["fname"],'w') as outfile:
                json.dump(dat,outfile)
        
    def write_line_data(self):
        '''
        Write line data (dictionairy) as .json
        '''
        import json
        for dat in self.line_data:
            json.dumps(dat,sort_keys=True)
            with open(dat["fname"],'w') as outfile:
                json.dump(dat,outfile)
