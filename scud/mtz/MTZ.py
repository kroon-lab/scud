#!/usr/bin/env cctbx.python

from __future__ import division
from iotbx import mtz
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
import numpy as np

class MTZClass(object):
#class MTZClass:

    def __init__(self,fname, array_name=None):
        self.name = fname
        self.array_name = array_name
        self.read_mtz_and_select_column()

    def hline_mtz(self,supercell_size=1):
        '''
        return list containing l=0 k=0 values
        '''
        diffArray = self.miller_array
        min_max = diffArray.min_max_indices()
        min_index, max_index  = min_max[0][0], min_max[1][0]+1
        selection_indices = [(x,0,0) for x in range(min_index,max_index) for y in range(min_index,max_index)]
        l_zero = diffArray.select_indices(selection_indices)
        dat = l_zero.data()
        result = []
        for i,miller in enumerate(l_zero.indices()):
            result.append([miller[0]/float(supercell_size),dat[i]])
        return result
        
    def make_P1(self):
        '''
        Convert to P1
        '''
        return self.miller_array.expand_to_p1()
        
    def read_mtz(self):
        '''
        read mtz file, convert to miller arrays and extract unit cell
        '''
        mtz_object = mtz.object(self.name)
        self.miller_array = mtz_object.as_miller_arrays()
        self.uc = self.miller_array[0].crystal_symmetry().unit_cell()
        self.symmetry = self.miller_array[0].crystal_symmetry()

    def read_mtz_and_select_column(self):
        '''
        Read mtz file, convert to object, when array_nun = None, print column name list
        and exit, else only load millers and chosen column.
        '''
        mtz_object = mtz.object(self.name)
        labs = mtz_object.column_labels()
        if self.array_name == None:
            print 'Choose column name:'
            print labs
            quit()
        else:
            ct = mtz_object.column_types()[labs.index(self.array_name)]
            column = mtz_object.get_column(self.array_name)
            column.set_type(ct)
            crystal = column.mtz_crystal()
            miller_set = crystal.miller_set(anomalous_flag=False)
            if ct == 'J':
                self.miller_array = miller_set.array(column.extract_values().as_double()).set_observation_type_xray_intensity()
            else:
                self.miller_array = miller_set.array(column.extract_values().as_double()).set_observation_type_xray_amplitude()
            self.symmetry = crystal.crystal_symmetry()
            self.uc = crystal.crystal_symmetry().unit_cell()
        
    def renumber_indices(self,supercell_size=None):
        """
        Renumber selected miller indices
        """
        # Create new space group and symmetry
        new_symmetry = crystal.symmetry(
            unit_cell = tuple(np.append(np.array(self.uc.parameters()[0:3])/supercell_size,
                                        np.array(self.uc.parameters())[3:7])),
            space_group_symbol="P1")
        # divide supercell miller indices by supercell size
        new_miller = flex.miller_index([(int(x/supercell_size),int(y/supercell_size),int(z/supercell_size))
                                        for x,y,z in miller_array.indices()])
        # create new miller set
        ms = miller.set(
            crystal_symmetry = new_symmetry,
            anomalous_flag = False,
            indices = new_miller)
        # create new miller array
        new_miller_array = miller_array.customized_copy(miller_set=ms,data=miller_array.data())
        return new_miller_array

    def select_indices(self,mtz_out_name=None,supercell_size=None,index_selection=None,make_P1=False):
        '''
        Selects indices according to supercell size
        '''
        if make_P1:
            newArray = self.make_P1()
        if index_selection == None:
            min_max = self.miller_array.min_max_indices()
            min_index, max_index  = min_max[0][0], min_max[1][0]+1
            hk_range = [i for i in range(0,max_index,supercell_size)] + [i for i in range(0,min_index,-supercell_size)]
            l_range = [i for i in range(0,max_index,supercell_size)]
            index_selection = [(h,k,l) for h in hk_range for k in hk_range for l in l_range]
        selection = self.miller_array.select_indices(index_selection)
        renumbered_selection = self.renumber_indices(miller_array=selection,
                                                     supercell_size = supercell_size)
        if mtz_out_name != None:
            self.write_array_as_mtz(renumbered_selection,
                                    mtz_out_name=mtz_out_name)
        return renumbered_selection,index_selection

    def show_miller_arrays(self):
        '''
        Print array names, and return choosen array number 
        '''
        return [ar.info().labels[0] for ar in self.miller_arrays]
    
    def slice_mtz(self,
                  slice_index='l'):
        '''
        Select an array from the MTZ file. From this select only the miller indices with l=0, write this as a MTZ file. Then apply centrosymmetry fill and numpy array, return this.

        Input:
        * self.mtz_object
        * mtz_out_name

        Returns:
        * numpy array with l=0 slice

'''
        # Miller array to slice:
        diffArray = self.miller_array
        # Find boundaries:
        min_max = diffArray.min_max_indices()
        # Selection indices based on slice_index value
        if slice_index == 'h':
            selection_indices = [(0,y,z) for y in range(min_max[0][1],min_max[1][1])
                                 for z in range(min_max[0][2],min_max[1][2])]
        if slice_index == 'k':
            selection_indices = [(x,0,z) for x in range(min_max[0][0],min_max[1][0])
                                 for z in range(min_max[0][2],min_max[1][2])]
        if slice_index == 'l':
            selection_indices = [(x,y,0) for x in range(min_max[0][0],min_max[1][0])
                                 for y in range(min_max[0][1],min_max[1][1])]
        # Select indices
        sel = diffArray.select_indices(selection_indices)
        # Create empty 2D array
        if slice_index == 'h':
            slice_2D = np.empty(((min_max[1][1]-min_max[0][1])+1,2*(min_max[1][2]-min_max[0][2])+1))
            slice_2D[:] = np.nan
            for sf in sel:
                # x,y coordinates in array
                x1 = sf[0][1] + min_max[1][1]
                y1 = sf[0][2] + min_max[1][2]
                # Assign intensity / structure factor to array index
                slice_2D[x1][y1] = sf[1]
                # Centro Symmetry:
                slice_2D[min_max[1][1] - ( x1 - min_max[1][1]),
                         min_max[1][2] - ( y1 - min_max[1][2])] = sf[1]
        if slice_index == 'k':
            slice_2D = np.empty(((min_max[1][0]-min_max[0][0])+1,2*(min_max[1][2]-min_max[0][2])+1))
            slice_2D[:] = np.nan
            for sf in sel:
                # x,y coordinates in array
                x1 = sf[0][0] + min_max[1][0]
                y1 = sf[0][2] + min_max[1][2]
                # Assign intensity / structure factor to array index
                slice_2D[x1][y1] = sf[1]
                # Centro Symmetry
                slice_2D[min_max[1][0] - ( x1 - min_max[1][0]),
                         min_max[1][2] - ( y1 - min_max[1][2])] = sf[1]
        if slice_index == 'l':
            slice_2D = np.empty(((min_max[1][0]-min_max[0][0])+1,(min_max[1][1]-min_max[0][1])+1))
            slice_2D[:] = np.nan
            for sf in sel:
                # x,y coordinates in array
                x1 = sf[0][0] + min_max[1][0]
                y1 = sf[0][1] + min_max[1][1]
                # Assign intensity / structure factor to array index
                slice_2D[x1][y1] = sf[1]
                # Centro Symmetry:
                slice_2D[min_max[1][0] - ( x1 - min_max[1][0]),
                min_max[1][1] - ( y1 - min_max[1][1])] = sf[1]
        # Return array
        return slice_2D

    def write_array_as_mtz(self,array=None,
                           label='SLC',
                           mtz_out_name=None):
        '''
        Write miller array as MTZ file
        '''
        if mtz_out_name == None: mtz_out_name = 'new.mtz'
        mtz_dataset = array.as_mtz_dataset(column_root_label=label,column_types='J')
        mtz_dataset.mtz_object().write(mtz_out_name)
 
def mtz_correlation(mtz_1 = None,
                    mtz_2 = None,
                    super_1 = 1,
                    super_2 = 1,
                    l = None):
    '''
    Calculate R-values and CC's between two mtz class objects
    '''

###########################################################################
#                   Expand Both sets to P1 if needed                      #
###########################################################################

    l.show_info('space group 1: {}'.format(mtz_1.symmetry.space_group_info()))
    l.show_info('space group 2: {}'.format(mtz_2.symmetry.space_group_info()))

###########################################################################
#                 Select indices and calculate R-Values                   #
###########################################################################

    l.process_message('Filtering and renumbering desired reflections...')
    # Check if spacegroup = P1, else expand to P1
    if str(mtz_1.symmetry.space_group_info()) != 'P 1':
        l.show_info('expanding mtz_1 to P1')
        mtz_1.make_P1()
    if str(mtz_2.symmetry.space_group_info()) != 'P 1':
        l.show_info('expanding mtz_2 to P1')
        mtz_2.make_P1()
    # If supercell size is not 1, renumber miller indices
    if super_1 != 1:
        mtz_1.renumber_indices(super_1)
    if super_2 != 1:        
        mtz_2.renumber_indices(super_2)
    # Get sets from both
    sel_1 = mtz_1.miller_array.set()
    sel_2 = mtz_2.miller_array.set()

    # Make sure both sets have the same reflections
    l.process_message('Selecting similar miller indices...')
    select_mtz_1 = mtz_1.miller_array.common_set(sel_1)
    select_mtz_2 = mtz_2.miller_array.common_set(sel_1)
    print select_mtz_1.size()
    print select_mtz_2.size()
    print select_mtz_1.data().as_numpy_array()[0:10]
    print select_mtz_2.data().as_numpy_array()[0:10]
    # Convert intensities to amplitudes
    l.process_message('Converting to complex amplitude array')
    select_mtz_1_amp = select_mtz_1.as_amplitude_array().set_observation_type_xray_amplitude()
    select_mtz_2_amp = select_mtz_2.as_amplitude_array().set_observation_type_xray_amplitude()
    tst1 = select_mtz_1_amp.data().as_numpy_array()
    tst2 = select_mtz_2_amp.data().as_numpy_array()
    print tst1[:10]
    print tst2[:10]
    tst_sf =  np.sum(tst1) / np.sum(tst2)
    print np.sum((tst1==0))
    print np.sum((tst2==0))
    
    print tst_sf
    l.process_message('Calculating R value...')
    # Scale factor calculation
    scale_factor = select_mtz_1_amp.scale_factor(select_mtz_2_amp.randomize_phases())
    # Calculate R-value
    print scale_factor
    r_value = select_mtz_1_amp.r1_factor(select_mtz_2_amp,
                                         scale_factor=scale_factor)
    r_value_t = select_mtz_1_amp.r1_factor(select_mtz_2_amp,
                                           scale_factor=tst_sf)
#                                         scale_factor=1)
    l.show_info('r_value: {:.2f}'.format(r_value))
    print 'tim',r_value_t

###########################################################################
#                            Calculate CC values                          #
###########################################################################

    l.process_message('Calculating overall CC value...')
    overall_CC = select_mtz_1_amp.correlation(select_mtz_2_amp).coefficient()
    l.show_info('Overall correlation: {:.2f}'.format(overall_CC))
    l.process_message('Calculating CC value per resolution shell...')
    # Set up binners
    select_mtz_1_amp.setup_binner(n_bins=20)
    select_mtz_2_amp.use_binning_of(select_mtz_1_amp)
    # Loop over resolution bins and calculate CC values
    for i_bin in select_mtz_1_amp.binner().range_all():
        sel_1 = select_mtz_1_amp.binner().selection(i_bin)
        sel_2 = select_mtz_2_amp.binner().selection(i_bin)
        sel_1_bin = select_mtz_1_amp.select(sel_1)
        sel_2_bin = select_mtz_2_amp.select(sel_2)
        cc_bin = sel_1_bin.correlation(sel_2_bin).coefficient()
        # Determine resolution range for bin
        legend = select_mtz_1_amp.binner().bin_legend(i_bin,show_counts=False)
        l.show_info('{} {}'.format(legend,cc_bin))
