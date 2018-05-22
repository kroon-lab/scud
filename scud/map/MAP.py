#!/usr/bin/env cctbx.python

from __future__ import division

import numpy as np
import os
from random import randint

import iotbx.ccp4_map
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import maptbx
from cctbx.array_family import flex
from scud.plt.Plot import Plot

class MapMapClass(object):
    '''
    Class used to compare maps with each other
    higher in heirachy then MapClass 

    input 2 maps
    '''
    def __init__(self,
                 map_1 = None,
                 map_1_dat = False,
                 map_2 = None,
                 map_2_dat = False,
                 p = None,
                 l = None):

        # Fisrt map (MAPClass)
        self.map_1 = map_1
        # Is measured
        self.map_1_dat = map_1_dat
        # First map (MAPClass)
        self.map_2 = map_2
        # Is measured
        self.map_2_dat = map_2_dat
        # Create prefix for map2map calculations
        self.prefix_generator(p = p,
                              l = l)
        # Create and / or check output folder
        if not os.path.isdir(self.prefix):
            os.mkdir(self.prefix)
        self.prefix = self.prefix+'/'
        # Data array map 1
        self.dat_1 = None
        # Data array map 2
        self.dat_2 = None
        # Data array difference (after scaling)
        self.diff = None
        # Template to write data arrays into
        self.template = None
        # CC-value between maps
        self.CC = None

    def bragg_mask(self,
                   p = None,
                   l = None):
        '''
        Returns a mask object that masks averything except where
        the original Bragg reflections have been / are.
        depends on size_h, size_k, size_l

        THIS WORKS NOW: add for loop for readability
        '''
        l.process_message('Creating Bragg mask')
        # Super cell sizes:
        s_h,s_k,s_l = p.params.size_h, p.params.size_k, p.params.size_l

        # grid map:
        grd = self.template.grd
        
        # Create all list with all Bragg positions as index (shifted)
        half_grd = (np.array(grd) - 1) / 2.
        
        # Future boolean arrays
        x_m,y_m,z_m = np.zeros((grd),dtype=bool),np.zeros((grd),dtype=bool),np.zeros((grd),dtype=bool)

        # # Creating boolean 1d list
        # Offset is for in between bragg locations
        n_h = np.arange(0,grd[0])
        x = ((n_h - half_grd[0]) % s_h == 0)
        n_k = np.arange(0,grd[1])
        y = ((n_k - half_grd[1]) % s_k == 0)
        n_l = np.arange(0,grd[2])
        z = ((n_l - half_grd[2]) % s_l == 0)

        # Setting 'lines' to be false
        x_m[x,:,:] = True
        y_m[:,y,:] = True
        z_m[:,:,z] = True

        # Create Bragg grid
        b_bragg = (x_m * y_m * z_m)

        # Convert to flex bool mask
        self.not_bragg_mask = flex.bool(True - b_bragg)
        
    def calc_CC(self,
                p = None,
                l = None):
        '''
        Calculate cross correlations between two maps

        THIS WORKS NOW: add for loop in bragg_CC partfor readability
        '''
        l.process_message('Calculating CC value...')
        # Converting (already masked) data to flex 1d doubles
        m1 = self.dat_1.as_double().as_1d()
        m2 = self.dat_2.as_double().as_1d()
        # Select values based on measured intensities or not (boolean mask)
        vals = ( m1 != -1000 )
        # Calculate liniear correlation between filtered maps
        mc = flex.linear_correlation(x = m1.select(vals),
                                     y = m2.select(vals))
        # Extract correlation coefficient from linear correlation object
        self.CC = mc.coefficient()
        # Calculate CC on Bragg and inter Bragg loactions (filter out rest)
        if p.params.bragg_cc:

            # Write Bragg and inter-Bragg data to file
            write_masked = False
            
            # Create Bragg related masks
            self.bragg_mask(p = p,
                            l = l)
            
            # Set all masked values to -1000 and make 1d double flex array of data
            grd = self.template.grd
            m_1_brg = self.dat_1.deep_copy().as_1d().set_selected(self.not_bragg_mask.as_1d(),-1000)
            m_1_brg.reshape(flex.grid(grd))
#            m_1_brg = self.dat_1.deep_copy().set_selected(self.not_bragg_mask,p.params.ignore_exp)
            if write_masked:
                self.template.write_map(dat = m_1_brg,
                                        name = self.prefix+'bragg_m_1_brg.map',
                                        l = l)
            m_1_brg = m_1_brg.as_double().as_1d()
#            m_2_brg = self.dat_2.deep_copy().set_selected(self.not_bragg_mask,p.params.ignore_exp)
            m_2_brg = self.dat_2.deep_copy().as_1d().set_selected(self.not_bragg_mask.as_1d(),-1000)
            m_2_brg.reshape(flex.grid(grd))
            if write_masked:
                self.template.write_map(dat = m_2_brg,
                                        name = self.prefix+'bragg_m_2_brg.map',
                                        l = l)
            m_2_brg = m_2_brg.as_double().as_1d()
            # Create 1d boolean of every location that is masked or not
            m_x_brg_vals = ( m_1_brg != p.params.ignore_exp )
            # Calculate actual correlation coefficents
            mc_brg = flex.linear_correlation(x = m_1_brg.select(m_x_brg_vals),
                                             y = m_2_brg.select(m_x_brg_vals))
            
            # Exctact and assign Correlation Coefficients:
            self.brg_CC = mc_brg.coefficient()
        # Print and save CC values
        l.show_info('overal CC-value     :  {}'.format(self.CC))
        if p.params.bragg_cc:
            l.show_info('Bragg CC-value      :  {}'.format(self.brg_CC))
        with open(self.prefix+'CC.txt','w') as f:
            print >> f, self.CC
        if p.params.bragg_cc:
            with open(self.prefix+'CC_brg.txt','w') as f:
                print >> f, self.brg_CC
            
    def create_mask(self,
                    p = None,
                    l = None):
        '''
        Create mask based on padding values of both input maps 
        If experimental masking value is different then when calculate
        '''
        l.process_message('Creating mask based on padding values...')
        # Create masks
        if p.params.map_1_exp:
            self.mask = (self.dat_1 == p.params.ignore_exp).set_selected( self.dat_2 == 0.0, True)
        if p.params.map_2_exp:
            self.mask = (self.dat_2 == p.params.ignore_exp).set_selected(self.dat_1 == 0.0, True)
        if p.params.map_1_exp == False and p.params.map_2_exp == False:
            self.mask = (self.dat_1 == 0.0).set_selected(self.dat_2 == 0.0, True)
        if p.params.map_1_exp == True and p.params.map_2_exp == True:
            self.mask = (self.dat_1 == p.params.ignore_exp).set_selected(self.dat_2 == p.params.ignore_exp, True)
        # Mask arrays
        self.dat_1.set_selected(self.mask,p.params.ignore_exp)
        self.dat_2.set_selected(self.mask,p.params.ignore_exp)
            
    def diff_map(self,
                 p = None,
                 l = None):
        '''
        Calculate Difference map between maps and write out
        '''
        l.process_message('Creating difference map...')
        self.diff = self.dat_1 - self.dat_2

        self.histo(p = p,
                   l = l,
                   dat = self.diff,
                   fname = self.prefix+'histo_difference_map.eps')
        if p.params.make_diff_maps:
            # Output map name
            diff_map_name = self.prefix+'difference.map'
            # Write diff map
            l.show_info('Writing difference map: {}'.format(diff_map_name))
            self.template.write_map(dat = self.diff,
                                    name = diff_map_name,
                                    l = l)

    def histo(self,
              p = None,
              l = None,
              dat = None,
              fname = None):
        '''
        Create histogram of input data array (random values)
        '''
        pl = Plot()
        # Check value distribution of difference map
        l.process_message('generating histogram...')
        limits = dat.all()
        m = 200000
        ind_list = zip(np.random.randint(limits[0]-1,size=m),
                       np.random.randint(limits[1]-1,size=m),
                       np.random.randint(limits[2]-1,size=m))
        l_1 = []
        for inds in ind_list:
            if not self.mask[inds]:
                l_1.append(dat[inds])
        # Create histogram
        pl.histogram(outname=fname,
                     dat = l_1,
                     bin_num= 50)

    def linear_CC(self,
                  p = None,
                  l = None,
                  m1 = None,
                  m2 = None,
                  vals = None):
        '''
        Wrapper for flex linear correlation coefficient function
        '''
        # Calculate liniear correlation between filtered maps
        mc = flex.linear_correlation(x = m1.select(vals),
                                     y = m2.select(vals))
        # Extract correlation coefficient from linear correlation object
        return mc.coefficient()
        
    def prefix_generator(self,
                         p = None,
                         l = None):
        '''
        Create prefix for output files
        '''
        pre_prefix = 'map2map_'
        if '/' in p.input.map_1: 
            head, sep, p_1 = p.input.map_1.partition('/')
        else:
            p_1 = p.input.map_1
        if '/' in p.input.map_2: 
            head, sep, p_2 = p.input.map_2.partition('/')
        else:
            p_2 = p.input.map_2
        self.prefix = pre_prefix+p_1.replace('.map','')+'_'+p_2.replace('.map','')

    def project_diff(self,
                     mask_value = -1000,
                     p = None,
                     l = None):
        '''
        Create a projection along l from the whole difference map
        Plot with 3 color color map where middle color is zero

        !!!!! Copy paste of map project function, should place somewhere else!
        '''
        l.process_message('Projecting difference map')
        # get data array
        dat = self.diff
        # Shape of array
        shape = dat.all()
        # Initiate 2D array of shape h,k
        projection = np.zeros(shape[0:2])
        # Counter array
        cnt = np.zeros(shape[0:2])
        # Convert array to numpy array
        dat_np = dat.as_numpy_array()
        l.show_info('Map converted to numpy array, projecting...')
        # 
        for i in xrange(shape[2]):
            slc = dat_np[:,:,i]
            cnt += (slc!=mask_value).astype(int)
            slc[slc==mask_value] = 0
            projection += slc
        cnt[cnt==0] = np.inf
        projection = projection / cnt #shape[2]
        # Projection statistics:
        l.process_message('Projection statistics:')
        l.show_info('Projection minimum: {}'.format(np.min(projection)))
        l.show_info('Projection maximum: {}'.format(np.max(projection)))
        l.show_info('Projection mean:    {}'.format(np.mean(projection)))
        pl = Plot()
        pl.projection(dat = projection,
                      outname = self.prefix+'projection_difference_map.eps',
                      diff = True,
                      l = l)
    def radial(self,
               p = None,
               l = None):
        '''
        Create and write out binned radial average plot of both maps
        '''
        l.process_message('Constructing radial profiles...')
        # Prep bin sizes based on resolution and number of bins
        d = self.template.reso[0] # Max resolution
        l.show_info('Max resolution: {:.2f}'.format(d))
        n = p.params.radial_bins # Number of bins
        # resolution list to bin in (1/d^3):
        d_pwr = np.power(list(np.arange(0,(1/d**3),(1/d**3)/(1.0*n))) + [1/d**3],(1/3.))

        x_axis = (d_pwr[1:] + d_pwr[:-1])/2. # X-axis used for plotting (later)

        # Different formats of grid needed, 1d, 3d and flex.grid
        grd = self.template.grd
        grd_1d = np.product(grd)
        grd_reshape = flex.grid(grd)

        # Needed for gridding around sites:
        unit_cell = self.template.unit_cell
        uc = self.template.uc
        # Not measured values
        invalid_mask = self.mask

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
        radial_map_1 = self.dat_1.deep_copy()
        radial_map_2 = self.dat_2.deep_copy()
        
        # Will become aniso maps:
        m1 = self.dat_1
        m2 = self.dat_2

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
            m2_bin_I = m2.select(total_mask.as_1d())
            # Calculate mean I
            m1_I_mean = flex.mean(m1_bin_I)
            m2_I_mean = flex.mean(m2_bin_I)

            # Calculate CC between maps for resolution sheel
            reso_CC = flex.linear_correlation(x = m1_bin_I.as_double(),
                                              y = m2_bin_I.as_double()).coefficient()

            # Append for later plotting
            reso_CC_ls.append([x_axis[i],reso_CC,(d_reso-d_star_min)])

            # Print reso table
            stat_string = '{:.2f}-{:.2f} ({} / {}) : {:.2f} / {:.2f}   CC: {:.2f}'.format(d_star_min,d_reso,m1_bin_I.size(),m2_bin_I.size(),m1_I_mean,m2_I_mean,reso_CC)
            stat_string_ls.append(stat_string)
            l.show_info(stat_string) # Write to file at some point!
            # Append mean I for shell
            m1_I_l.append(m1_I_mean)
            m2_I_l.append(m2_I_mean)
            
            # Subtract radial average from masked part of aniso_map:
            radial_map_1.set_selected(total_mask,m1_I_mean)
            radial_map_2.set_selected(total_mask,m2_I_mean)
            
            # Move upper mask to lower mask of next shell
            mask_lower_binary = mask_upper_binary
            # Update lower resolution limit
            d_star_min = d_reso

        # Append to I_list container
        I_lists.append(m1_I_l)
        I_lists.append(m2_I_l)

        # Subtract radial mask from dat set invalid to 0 and write map
        l.process_message('Creating and writing aniso map...')
        
        # Create aniso maps
        #if not p.params.map_1_exp:
        m1 = m1 - radial_map_1
        m2 = m2 - radial_map_2
        
        # Final aniso maps invalid mask set to 0 because intensities are distributed around 0:
        m1.set_selected(invalid_mask,0)
        m2.set_selected(invalid_mask,0)
        
        m_name, radial_ls = [p.input.map_1,p.input.map_2], [radial_map_1,radial_map_2]

        for m,mp in enumerate([m1,m2]):
            self.template.write_map(dat = mp,
                                    name = self.prefix+m_name[m][:-4]+'_aniso.map',
                                    l = l)
            # Double check, write the map containing radial averages as wel.
            self.template.write_map(dat = radial_ls[m],
                                    name = self.prefix+m_name[m][:-4]+'_radial.map',
                                    l = l)
        
        # Calculate CC_aniso
        l.process_message('Calculating aniso CC')
        # prep maps, set invalid to -1000, mask double and 1d
        m1 = m1.set_selected(invalid_mask,-1000).as_double().as_1d()
        m2 = m2.set_selected(invalid_mask,-1000).as_double().as_1d()
        # Calc indices used for CC calculation
        vals = (m1 != -1000)        
        # Calc linear CC
        aniso_CC = self.linear_CC(p = p,l = l,
                                  m1 = m1,
                                  m2 = m2,
                                  vals = vals)
        l.show_info('Aniso CC value: {}'.format(aniso_CC))

        # Write aniso_CC to file
        with open(self.prefix+'CC_aniso.txt','w') as f:
                print >> f, aniso_CC

        # Write reso table to file:
        with open(self.prefix+'radial_table.txt','w') as f:
            for s in stat_string_ls: print >> f, s

        # IMPLEMENT 'FLAT' / 'SHAPE' CC
        m1_aniso_pos = ( m1 > 0 ).set_selected( m1 == -1000, False)
        m1_aniso_neg = ( m1 < 0 ).set_selected( m1 == -1000, False)
        m2_aniso_pos = ( m2 > 0 ).set_selected( m2 == -1000, False)
        m2_aniso_neg = ( m2 < 0 ).set_selected( m2 == -1000, False)

        m1.set_selected(m1_aniso_neg,-1).set_selected(m1_aniso_pos,1)
        m2.set_selected(m2_aniso_neg,-1).set_selected(m2_aniso_pos,1)
        
        vals = ( m1 != -1000 )
        
        shape_CC = self.linear_CC(p = p,l = l,
                                  m1 = m1,
                                  m2 = m2,
                                  vals = vals)
        l.show_info('Shape CC value: {}'.format(shape_CC))

        m1.reshape(grd_reshape)
        m1 = m1.set_selected(invalid_mask, 0)
        m2.reshape(grd_reshape)
        m2 = m2.set_selected(invalid_mask, 0)
        for m,mp in enumerate([m1,m2]):
            self.template.write_map(dat = mp   ,
                                    name = self.prefix+m_name[m][:-4]+'_shape.map',
                                    l = l)        
        
        # Plot CC-reso values in bar graph
        x,y,w = map(list, zip(*reso_CC_ls))
        pt = Plot()
        pt.collect_bar_data(xlabel='Reso ($\AA^{-1}$)',
                            ylabel='CC',
                            label='CC {} vs {}'.format(m_name[0][:-4],m_name[1][:-4]),
                            color=(0.2,0.2,0.8),
                            x=x,
                            y=y,
                            width=w,
                            fname=self.prefix+'reso_CC.eps')
        pt.write_bar_data()
        pt.bar_plot(plotfile_name=self.prefix+'reso_CC.eps')
        
        # Plotting Radial profiles
        l.process_message('Writing Line data and plotting...')
        f_1,f_2 = m_name[0][:-4], m_name[1][:-4]
        colors = [(0.2,0.2,0.8),(0.8,0.8,0.2)] # Plot colors (RGB)
        nm = [self.prefix+f_1+'_radial.json',
              self.prefix+f_2+'_radial.json'] # Plot names
        ln = [f_1,f_2] # Line labels
        pt = Plot() # call Plot class
        for i,I_l in enumerate(I_lists):
            # Create line data for saving and plotting
            pt.collect_line_data(xlabel='Reso (1/d$^{3}$)',
                                 ylabel='Intensity',
                                 label=ln[i],
                                 color=colors[i],
                                 x=x_axis,
                                 y=I_l,
                                 fname=nm[i])
        # Write lines as .json
        pt.write_line_data()
        # Create final plot and write out
        pt.lines_plot(plotfile_name=self.prefix+'radial_profiles.eps')
        
    def resize_maps(self,
                    p = None,
                    l = None):
        '''
        Get two maps, get size difference, slice the bigger one
        '''
        l.process_message('Comparing map sizes...')
        l.show_info('map 1 grid: {}'.format(self.map_1.grd))
        l.show_info('map 2 grid: {}'.format(self.map_2.grd))
        l.show_info('\n')
        
        # Difference in map sizes and borders
        diffs = np.abs(np.array(self.map_1.grd) - np.array(self.map_2.grd)) 
        side_diffs = diffs/2
        l.show_info('grd difference: {}'.format(diffs))
        l.show_info('slice borders: {}'.format(side_diffs))

        l.process_message('Slicing bigger map')
        # if map_1 is bigger:
        if self.map_1.grd > self.map_2.grd:
            l.show_info('Slicing map {}'.format(p.input.map_1))
            # slice map
            self.dat_1 = self.map_1.maps.data[int(side_diffs[0]):self.map_1.grd[0]-int(side_diffs[0]),
                                              int(side_diffs[1]):self.map_1.grd[1]-int(side_diffs[1]),
                                              int(side_diffs[2]):self.map_1.grd[2]-int(side_diffs[2])]
            self.dat_2 = self.map_2.maps.data
            self.template = self.map_2
        if self.map_1.grd < self.map_2.grd:
            l.show_info('Slicing map {}'.format(p.input.map_2))
            # slice map
            self.dat_2 = self.map_2.maps.data[int(side_diffs[0]):self.map_2.grd[0]-int(side_diffs[0]),
                                              int(side_diffs[1]):self.map_2.grd[1]-int(side_diffs[1]),
                                              int(side_diffs[2]):self.map_2.grd[2]-int(side_diffs[2])]
            self.dat_1 = self.map_1.maps.data
            self.template = self.map_2
            # write maps (second map is written to overlay origins)
        if self.map_1.grd == self.map_2.grd:
            self.dat_1 = self.map_1.maps.data
            self.dat_2 = self.map_2.maps.data
            self.template = self.map_1
            l.process_message('Both maps have the same size... Continuing')
        self.dat_1.reshape(flex.grid(self.template.grd))
        self.dat_2.reshape(flex.grid(self.template.grd))

    def scale_maps(self,
                   p = None,
                   l = None,
                   lsq = True):
        '''
        Scale maps based on middle 50% of values in a QQ plot
        '''
        if not lsq:
            l.process_message('QQ-plot will be used to scale maps...')
            # Plot histograms of unscaled maps
            self.histo(p = p,
                       l = l,
                       dat = self.dat_1,
                       fname = self.prefix+'histo_map_1_unscaled.eps')
            self.histo(p = p,
                       l = l,
                       dat = self.dat_2,
                       fname = self.prefix+'histo_map_2.eps')
            # Scale maps
            pl = Plot()
            if p.params.scale_maps:
                l.process_message('Scaling Maps...')
                l.process_message('Creating QQ-plot...')
                limits = self.dat_1.all()
                ind_list = [[randint(0,limits[0]-1),
                             randint(0,limits[1]-1),
                             randint(0,limits[2]-1)]
                            for i in range(20000)]
                l.process_message('Random indices generated...')
                l_1, l_2 = [],[]
                for inds in ind_list:
                    if not self.mask[inds]:
                        l_1.append(self.dat_1[inds])
                        l_2.append(self.dat_2[inds])
                        # Check if values in list are measured/calculated values
                        # Sort lists
                l_1 = np.sort(np.array(l_1))
                l_2 = np.sort(np.array(l_2))
                l.show_info('{} values will be used for QQ-plot'.format(len(l_1)))
                l.process_message('Creating QQ-plot...')
                # Fit and plot
                xlab,ylab = '$I_{calc}$','$I_{calc}$'
                if p.params.map_1_exp:
                    xlab = '$I_{obs}$'
                if p.params.map_2_exp:
                    ylab = '$I_{obs}$'
                intersect, slope = pl.QQ_plot(dat_1=l_1,
                                              dat_2=l_2,
                                              xlab=xlab,
                                              ylab=ylab,
                                              out_name=self.prefix+'QQ_plot_out.eps',
                                              l = l)

        ###########################################################################
        #                 Actual scaling and writing new map                      #
        ###########################################################################

                if p.params.map_1_exp:
                    self.dat_1 = (self.dat_1*slope)+intersect
                    self.dat_1.set_selected(self.mask,p.params.ignore_exp)
                    m_1_name = self.prefix+'scaled.map'
                    l.show_info('Writing map: {}'.format(m_1_name))
                    self.template.write_map(dat = self.dat_1,
                                            name = m_1_name,
                                            l = l)
                    m_2_name = self.prefix+'compare.map'
                    l.show_info('Writing map: {}'.format(m_2_name))
                    self.template.write_map(dat = self.dat_2,
                                            name = m_2_name,
                                            l = l)
                if p.params.map_2_exp:
                    self.dat_2 = (self.dat_2*slope)+intersect
                    self.dat_2.set_selected(self.mask,p.params.ignore_exp)
                    m_2_name = self.prefix+'scaled.map'
                    l.show_info('Writing map: {}'.format(m_2_name))
                    self.template.write_map(dat = self.dat_2,
                                            name = m_2_name,
                                            l = l)
                    m_1_name = self.prefix+'compare.map'
                    l.show_info('Writing map: {}'.format(m_1_name))
                    self.template.write_map(dat = self.dat_1,
                                            name = m_1_name,
                                            l = l)
                    self.histo(p = p,
                               l = l,
                               dat = self.dat_1,
                               fname = self.prefix+'histo_map_1_scaled.eps')

        if lsq:
            l.process_message('Least squares will be used to scale maps:')
            # Prepare for least squares, map 1 will be scaled to map 2

            self.histo(p = p,
                       l = l,
                       dat = self.dat_1,
                       fname = self.prefix+'histo_lsq_map_1_unscaled.eps')
            
            np_mask = True - self.mask.as_numpy_array().flatten()
            
            np_dat_1 = self.dat_1.as_numpy_array().flatten()
            np_vals_1 = np_dat_1[np_mask]
            
            np_dat_2 = self.dat_2.as_numpy_array().flatten()
            np_vals_2 = np_dat_2[np_mask]

            vals_1_stack = np.vstack([np_vals_1,np.ones(len(np_vals_1))]).T

            # use numpy linear least squares function to determine multiplier and constant
            l.process_message('Performing LSQ...')
            # Polyfit
            m,c = np.linalg.lstsq(vals_1_stack,np_vals_2)[0]
            l.show_info('LSQ, y=ax+b, Multiplier(a)={}, constant(b)={}'.format(m,c))

            # Implement somewhere else (important plot):
            # import matplotlib.pyplot as plt
#             print 'Scatter for nick'
            
#             idx = np.random.randint(0,len(np_vals_1),100000)
#             plt.scatter(np_vals_1[idx] ,np_vals_2[idx])
#             x = np_vals_1[idx]
#             y = m*x + c
#             plt.plot(x,y,'-')
#             plt.show()
# #            quit()
            
            # Scale map 1
            self.dat_1 = ( self.dat_1 * m ) + c
            # Write scaling results
            self.dat_1.set_selected(self.mask,p.params.ignore_exp)
            m_1_name = self.prefix+'lsq_scaled.map'
            l.show_info('Writing map: {}'.format(m_1_name))
            self.template.write_map(dat = self.dat_1,
                                    name = m_1_name,
                                    l = l)
            m_2_name = self.prefix+'lsq_compare.map'
            l.show_info('Writing map: {}'.format(m_2_name))
            self.template.write_map(dat = self.dat_2,
                                    name = m_2_name,
                                    l = l)
            # Create histogram from scaled map
            self.histo(p = p,
                       l = l,
                       dat = self.dat_1,
                       fname = self.prefix+'histo_lsq_map_1_scaled.eps')
            self.histo(p = p,
                       l = l,
                       dat = self.dat_2,
                       fname = self.prefix+'histo_lsq_map_2.eps')
            
        
class MAPClass(object):

    def __init__(self,fname):
        self.name = fname
        self.read_map()

    def ind_to_d_2(self,ind1,ind2):
        '''
        Returns the resolution in A^-2 of a map index
        '''
        # Index (origin corrected) to A^-2
        recip_coord = (ind2-self.origin_shift)*self.voxs
        # return distance in A^-2
        return np.linalg.norm(ind1-recip_coord)**2

    def map_info(self,l=None):
        '''
        Print map info to screen and log, 
        calculate/save extra values needed later
        '''

###########################################################################
#                            Extra Values                                 #
###########################################################################
        
        self.grd = self.maps.unit_cell_grid
        self.uc = self.maps.unit_cell_parameters
        self.unit_cell = self.maps.unit_cell()
        self.reso = 1 / (np.array(self.uc[0:3])/2)
        self.voxs = np.array(self.uc[0:3])/np.array(self.grd)
        self.origin_shift = ((np.array(self.grd) - 1) / 2) + 1
        self.voxel_size = 1/self.voxs
        
###########################################################################
#                            Show Info                                    #
###########################################################################
        
        l.show_info('header_min: {}'.format(self.maps.header_min))
        l.show_info('header_max: {}'.format(self.maps.header_max))
        l.show_info('header_mean: {}'.format(self.maps.header_mean))
        l.show_info('header_rms: {}\n'.format(self.maps.header_rms))
        l.show_info('origin: {}'.format(self.maps.data.origin()))
        l.show_info('origin shift: {}'.format(self.origin_shift))
        l.show_info('grid: {}'.format(self.grd))
        l.show_info('abc (A^-1): {}'.format(self.uc[0:3]))
        l.show_info('angles (deg): {}'.format(self.uc[3:7]))
        l.show_info('Space Group Number: {}'.format(self.maps.space_group_number))
        l.show_info('Resolution (A): {}'.format(self.reso))
        l.show_info('voxel size (A^-1): {}'.format(self.voxs))
        l.show_info('voxel size (A): {}'.format(1/self.voxs))
    
    def max_reso_d2(self):
        '''
        Return max resolution (A^-2) of map

        map is a cube, while reciprocal space is a sphere so
        cornes cube do not count
        '''
        return self.ind_to_d_2(np.array([0,0,0]),np.array([self.grd[0],
                                                           self.origin_shift[1],
                                                           self.origin_shift[2]]))

    def projection(self,
                   l = None,
                   outname = 'project.eps',
                   mask_value = -1000):
        '''
        Create a projection along l from the whole map
        '''
        l.process_message('Projecting difference map')
        # get data array
        dat = self.maps.data
        # Shape of array
        shape = dat.all()
        # Initiate 2D array of shape h,k
        projection = np.zeros(shape[0:2])
        # Counter array
        cnt = np.zeros(shape[0:2])
        # Convert array to numpy array
        dat_np = dat.as_numpy_array()
        l.show_info('Map converted to numpy array, projecting...')
        # 
        for i in xrange(shape[2]):
            slc = dat_np[:,:,i]
            cnt += (slc!=mask_value).astype(int)
            slc[slc==mask_value] = 0
            projection += slc
        cnt[cnt==0] = np.inf
        projection = projection / cnt #shape[2]
        # Projection statistics:
        l.process_message('Projection statistics:')
        l.show_info('Projection minimum: {}'.format(np.min(projection)))
        l.show_info('Projection maximum: {}'.format(np.max(projection)))
        l.show_info('Projection mean:    {}'.format(np.mean(projection)))
        projection_txt = 'projection.txt'
        l.process_message('Writing projection to file: {}'.format(projection_txt))
        np.savetxt(projection_txt,projection)
        pl = Plot()
        pl.projection(dat = projection,
                      outname = outname,
                      diff = False,
                      l = l)

    def read_map(self):
        '''
        Read map
        '''
        self.maps = iotbx.ccp4_map.map_reader(file_name=self.name)

    def shift_origin(self,
                     l = None):
        '''
        Shift origin of map to 0,0,0
        '''
        new_name = 'shifted_map.map'
        map_data = self.maps.data.shift_origin()
        l.process_message('Shifting map origin')
        self.write_map(dat = map_data,
                       name = new_name,
                       l = l)
        l.process_message('Done')
        return new_name
        
    def subtract_radial(self,
                        l = None,
                        radial = None):
        '''
        Subtract radial average from map data
        '''
        l.process_message('Subtracting radial average from map...')
        
        # center of reciprocal space (A^-2)
        center = np.array([0,0,0])
        
        # Extract radial profile info from dictionary
        d2_range_list = radial["d2_range"]
        i_list = radial["i_mean"]
        
        # map data
        dat = self.maps.data.deep_copy()

        # Loop over all indices of the map
        for x in range(0,self.grd[0]):
            if x % 50 == 0: l.show_info('x = {}'.format(x))
            for y in range(0,self.grd[1]):
                for z in range(0,self.grd[2]):
                    # calculate hkl like value
#                    hkl = np.array([x,y,z]) - self.origin_shift
                    # calculate resolution (A^-2)
                    d2 = self.ind_to_d_2(center,np.array([x,y,z]))
                    # Find where hkl lies on radial profile
                    for i,d2_range in enumerate(d2_range_list):
                        if d2_range[0] < d2 < d2_range[1]:
                            dat[x,y,z] = dat[x,y,z] - i_list[i]
                        
        # Write map to file
        assert dat is not None
        self.write_map(dat = dat,
                       name = 'radial_corrected.map',
                       l = l)

    def write_map(self,dat = None,
                  name = None,
                  l = None):
        '''
        Write map to file
        If data not given, use self.maps.data
        if name not given, use new_map.map
        '''
        # if no input revert to defaults
        if dat == None: dat = self.maps.data
        if name == None: name = 'new_map.map'

        # write map to file
        iotbx.ccp4_map.write_ccp4_map(        
            file_name= name,
            unit_cell = uctbx.unit_cell(self.uc),
            space_group= sgtbx.space_group_info("P1").group(),
            map_data = dat.as_double(),
            labels=flex.std_string(["iotbx.ccp4_map.tst"]))
        l.process_message('Map written...')
