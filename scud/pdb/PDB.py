#!/usr/bin/env cctbx.python

from __future__ import division
import iotbx
from iotbx import pdb
import numpy as np
from scitbx.array_family import flex

from scud.general.method_library import datetime_string
from scud.general.method_library import vecDistSquare


class PDBClass(object):
    '''
    TODO: 
    * docstrings :P
    * ensemble trimmer based on either:
        * Structure Factor variance (example: lowest 10% R-value structures with respect to the average structure?)
        * RMSD, X structures with lowest RMSD compared to average structure?
    '''
    
    def __init__(self,
                 fname=None,
                 selection='All',
                 symmetry=None):
        self.fname = fname
        self.selection = selection + ' and (altid " " or altid A)'
        self.symmetry = symmetry
        self.hierarchy = None
        self.center_of_mass = None
        self.new_hierarchy = None
        self.read_pdb()
        self.renum_and_count_models()
        self.model_num = len(self.hierarchy.models())
        self.model_list = None
        self.target_b = None
        self.mask = None
        self.masked_target_b = None
        self.target_B_CA = None
        
    def average_B(self,
                  outlier_rejection=True,
                  m_sigma = 2,
                  make_mask = False,
                  dat = None):
        '''
        Return average B-factor value of PDB selected atoms
        '''
        if outlier_rejection:
            # Discard B-factors larger then m_sigma SD from the median of list 
            if dat == None: dat = np.array([atom.b for atom in self.hierarchy.atoms()])
            # Change to median deviation
            if make_mask:
                self.mask = (abs(dat-np.median(dat)) < m_sigma * np.std(dat))
                self.masked_target_b = dat[self.mask]
            return np.mean(dat[abs(dat - np.median(dat)) < m_sigma * np.std(dat)])
        else:
            # Just return the average
            return np.mean(np.array([atom.b for atom in self.hierarchy.atoms()]))

    def B_factor_to_RMSF(self,
                         b = None):
        '''
        Convert B-factor to RMSF in A
        '''
        return np.sqrt(b / (8 * (np.pi**2) ))
        
    def create_random_model_list(self,
                                 l = None,
                                 ensemble_size = None):
        '''
        Return a list with models ordered randomly, length depends on input
        or the amount of models in the input ensemble

        (for rigid_body.py)
        '''
        if self.model_num > 1:
            l.warning('Ensemble is used as input!!')
            l.process_message('Generating model list')
            if self.model_num > ensemble_size:
                l.warning('More models in input ensemble then desired output ensemble, changing output ensemble size')
                ensemble_size = self.model_num
            # how many times do all models fit in desired ensemble_size?
            n= ensemble_size // self.model_num
            model_list = np.array(range(self.model_num)*n)
            # Fill up with random models to desired ensemble size
            rest =  np.random.randint(self.model_num,size=ensemble_size-len(model_list))
            # Final model list
            model_list = np.concatenate([model_list,rest]).astype(int)
        else:
            # Generate array containing 0's, for 1 model PDB files
            model_list = np.array([0]*ensemble_size)
        return model_list

    def create_operation_list(self,
                              trans_only = False,
                              rot_only = False,
                              ensemble_size = None,
                              trans_sigma = None,
                              rot_sigma = None,
                              l = None ):
        '''
        Combine rotation and translation list depending on parameters
        '''
        if trans_only:
            rotation_list = [(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0) for i in xrange(ensemble_size)]
        translation_vector_list = self.create_translation_vector_list(mu=0,
                                                                      translation_sigma=trans_sigma,
                                                                      ensemble_size=ensemble_size)
        if not trans_only:
            # Generate rotation matrices based on normal distribution
            rotation_list = self.create_rotation_matrix_list(rotation_sigma_deg=rot_sigma,
                                                             ensemble_size=ensemble_size) 
        return translation_vector_list, rotation_list
    
    def create_rotation_matrix_list(self,rotation_sigma_deg,ensemble_size):
        '''
        Create rotation matrix list. 
        Generate random angles from normal distribution.
        Then generate random vector
        Calculate rotation matrix around the vector, using angle.
        '''
        rotation_sigma = np.deg2rad(rotation_sigma_deg)
        angle_list = np.random.normal(0.0,rotation_sigma,ensemble_size)
        rotation_matrix_list = []
        for angle in angle_list:
            xy = np.random.randint(-100,high=100,size=2)/100.
            z = 1 - xy[0]**2 - xy[1]**2
            xyz = np.array([xy[0],xy[1],z])
            rotation_matrix_list.append(self.rotation_matrix(xyz,angle))
        return rotation_matrix_list

    def create_translation_vector_list(self,mu,translation_sigma,ensemble_size):
        '''
        Create normal distributed (3D) translation vectors
        '''
        return np.random.multivariate_normal((mu,mu,mu),[[translation_sigma**2,0,0],[0,translation_sigma**2,0],[0,0,translation_sigma**2]],ensemble_size)

        
    def calculate_center_of_mass(self):
        '''
        Return Center of Mass from a selection/hierarchy
        '''
        return self.hierarchy.atoms().extract_xyz().mean()

    def detached_model_copy(self,n=0):
        '''
        Return a detached model copy of the hierarchy
        '''
        return self.hierarchy.models()[n].detached_copy()
    
    def ensemble_to_B_factor(self,ensemble_hierarchy=None,
                             ca_only=True):
        '''
        Calculate RMSF (CA) from an ensemble against the average position and convert to B-Factors
        '''
        sel_cache = ensemble_hierarchy.atom_selection_cache()
        if ca_only:
            selection_cache = sel_cache.selection('name CA and (altid " " or altid A)')
        else:
            selection_cache = sel_cache.selection('(altid " " or altid A)')
        ensemble_hierarchy = ensemble_hierarchy.select(selection_cache)
        for c,chain in enumerate(ensemble_hierarchy.models()[0].chains()):
            chainID = chain.id
            res_num = [int(resid.resseq) for resid in ensemble_hierarchy.models()[0].chains()[c].residues()]
            cnt = 0 
            for i,model in enumerate(ensemble_hierarchy.models()):
                if i == 0:
                    xyz = model.chains()[c].atoms().extract_xyz()
                else:
                    xyz+=model.chains()[c].atoms().extract_xyz()
                cnt+=1
            av_xyz = xyz*(1./cnt)
            for i,model in enumerate(ensemble_hierarchy.models()):
                if i == 0:
                    xyz = model.chains()[c].atoms().extract_xyz()
                    d_square = vecDistSquare(xyz,av_xyz)
                else:
                    xyz = model.chains()[c].atoms().extract_xyz()
                    d_square += vecDistSquare(xyz,av_xyz)
            av_d_square = d_square*(1./cnt)
            B = (8./3) * (np.pi**2) * av_d_square
        if ca_only:
            return np.array(res_num),np.array(B)
        else:
            return range(len(av_xyz)),np.array(B)

    def extract_ca_b_factors(self,single_structure_hierarchy=None):
        '''
        Xtract CA B-Factors from a single model PDB file and return a list (will be JSON) wih results
'''
        sel_cache = single_structure_hierarchy.atom_selection_cache()
        selection_cache = sel_cache.selection('name CA and (altid " " or altid A)')
        single_structure_hierarchy = single_structure_hierarchy.select(selection_cache)
        ca_b_factors = []
        x_data = []
        y_data = []
        for atom in single_structure_hierarchy.atoms():
            x_data.append(int(atom.fetch_labels().resid()))
            y_data.append(atom.b)
        return np.array(x_data), np.array(y_data)

    def minimum_B_Factor(self):
        '''
        Returns minimum Bfactor value
'''
        b_min = 10000
        for atom in self.hierarchy.atoms():
            if atom.b < b_min: b_min = atom.b
        return b_min

    def make_rot_trans_ens(self,
                           trans_only = False,
                           rot_only = False,
                           trans_sigma = None,
                           rot_sigma = None,
                           center_of_mass = None,
                           l = None):
        self.start_or_attach_hierarchy()
        #### Actual rotation/translation
        translation_vector_list, rotation_list = self.create_operation_list(trans_only = trans_only,
                                                                            rot_only = rot_only,
                                                                            ensemble_size=len(self.model_list),
                                                                            trans_sigma = trans_sigma,
                                                                            rot_sigma = rot_sigma,
                                                                            l = l)
        for i,(translation_vector,rotation_matrix) in enumerate(zip(translation_vector_list,rotation_list)):
            detached_model = self.detached_model_copy(n=self.model_list[i])
            rotated_model = self.rotate_model(detached_model,
                                              rotation_matrix,
                                              center_of_mass)
            translated_model = self.translate_model(rotated_model,
                                                    translation_vector)
            self.start_or_attach_hierarchy(detached_model=translated_model)
        
    def maximum_B_Factor(self):
        '''
        Returns maximum Bfactor value
'''
        b_max = 0
        for atom in self.hierarchy.atoms():
            if atom.b > b_max: b_max = atom.b
        return b_max
    
    def read_pdb(self):
        '''
        Reads a PDB file using iotbx.pdb

        Input:
        * filename

        Returns:
        * symmetry (global parameter)
        * pdb hierarchy (with or without selection filter)
        
        '''
        pdb_in = iotbx.pdb.input(file_name=self.fname)
        self.symmetry=pdb_in.crystal_symmetry()
        raw_hierarchy = pdb_in.construct_hierarchy()
        sel_cache = raw_hierarchy.atom_selection_cache()
        selection_cache = sel_cache.selection(self.selection)
        self.hierarchy = raw_hierarchy.select(selection_cache)

    def rechain(self, model=None):
        """
        Give chains in model unique identifiers
        """
        l = iotbx.pdb.systematic_chain_ids()
        for n,chain in enumerate(model.chains()):
            chain.id = l[n]
            
    def renum_and_count_models(self,
                               hierarchy = None):
        """
        Logically renumber models in an ensemble hierarchy
        """
        if hierarchy == None: hierarchy=self.hierarchy
        for i,model in enumerate(hierarchy.models()):
            model.id=str(i)
        self.model_num = int(hierarchy.models()[-1].id)

    def RMSF_to_B_factor(self,
                         rmsf = None):
        '''
        Convert B-factor to RMSF in A
        '''
        return (rmsf**2)* (8 * (np.pi**2) )
    
    def rotation_matrix(self,xyz,a):
        '''
        Create rotation matrix
        calc cos and sin first
        Then fill list in cctbx format with rotation matrix elements and return
'''
        c = np.cos(a)
        s = np.sin(a)
        R = (  c + ((xyz[0]**2)*(1-c)),
               ((xyz[0]*xyz[1])*(1-c))-(xyz[2]*s),
               ((xyz[0]*xyz[2])*(1-c))+(xyz[1]*s),
               ((xyz[1]*xyz[0])*(1-c))+(xyz[2]*s),
               c + ((xyz[1]**2)*(1-c)),
               ((xyz[1]*xyz[2])*(1-c))-(xyz[0]*s),
               ((xyz[2]*xyz[0])*(1-c))-(xyz[1]*s),
               ((xyz[2]*xyz[1])*(1-c))+(xyz[0]*s),
               c + ((xyz[2]**2)*(1-c)) )
        return R
        
    def rotate_model(self,model,rotation_matrix,center_of_mass):
        '''
        Rotate model. First use center of mass to translate to origin, then apply rotation matrix. After rotation use center of mass to back translate to original position.
        '''
        model.atoms().set_xyz( ( rotation_matrix * (model.atoms().extract_xyz() - center_of_mass ) ) + center_of_mass )
        return model

    def set_occ(self, occ=1.00):
        '''
        Set all B-Factors in a hierarchy to zero
        '''
        for atom in self.hierarchy.atoms():
            atom.occ = 1.00
    
    def set_B_zero(self):
        '''
        Set all B-Factors in a hierarchy to zero
        '''
        for atom in self.hierarchy.atoms():
            atom.b = 0.0

    def simplex_target_func(self,
                            trans_sigma = None,
                            rot_sigma = None,
                            center_of_mass = None,
                            trans_only = False,
                            rot_only = False,
                            l = None):
        '''
        Target function for simplex minimizer,
        create ensemble and calculate rms between original
        B-factor profile and B-factor from rms.
        '''
        # Create ensemble
        self.make_rot_trans_ens(trans_only = trans_only,
                                rot_only = rot_only,
                                trans_sigma = trans_sigma,
                                rot_sigma = rot_sigma,
                                center_of_mass = center_of_mass,
                                l = l)
        # Get b-factors
        atom,b = self.ensemble_to_B_factor(ensemble_hierarchy=self.new_hierarchy,
                                           ca_only = False)
        # Return rmsd with target b-factors
        return np.sqrt(np.mean((self.masked_target_b - b[self.mask])**2))

    def simplex_target_func_ca(self,
                               trans_sigma = None,
                               rot_sigma = None,
                               center_of_mass = None,
                               rb_type = None,
                               l = None):
        '''
        Target function for simplex minimizer,
        create ensemble and calculate rms between original
        B-factor profile and B-factor from rms.
        '''
        # Create ensemble
        if rb_type == 'mix':
            trans_only, rot_only = False, False
        if rb_type == 'trans':
            trans_only, rot_only = True, False
        if rb_type == 'rot':
            trans_only, rot_only = False, True
            
        self.make_rot_trans_ens(trans_only = trans_only,
                                rot_only = rot_only,
                                trans_sigma = trans_sigma,
                                rot_sigma = rot_sigma,
                                center_of_mass = center_of_mass,
                                l = l)
        # Get b-factors
        atom,b = self.ensemble_to_B_factor(ensemble_hierarchy=self.new_hierarchy,
                                           ca_only = True)
        # Return rmsd with target b-factors
        return np.sqrt(np.mean((self.target_B_CA - b)**2))
    
    def start_or_attach_hierarchy(self,
                                  detached_model=None,
                                  return_hierarchy=False):
        '''
        Either start an empty root hierarchy or append a model to it.
        '''
        if return_hierarchy: return self.new_hierarchy
        if detached_model == None: self.new_hierarchy = iotbx.pdb.hierarchy.root()
        else:
            self.new_hierarchy.append_model(detached_model)
    
    def translate_model(self,
                        input_detached_model,
                        translation_vector):
        '''
        Translate model and return it.
        '''
        input_detached_model.atoms().set_xyz( input_detached_model.atoms().extract_xyz() + (translation_vector) )
        return input_detached_model

    def trim_ensemble(self,
                      target = 100,
                      out_name = None,
                      skip_method = True,
                      l = None):
        '''
        Reduce the number of model in an ensemble
        '''
        l.process_message('Trimming ensemble...')
        if out_name == None: out_name='trimmed_ensemble_out.pdb'
        if skip_method:
            l.process_message('Using simple method for trimming')
            skip_num = int(self.model_num / target)
            l.show_info('skipping {} models'.format(skip_num))
            i,cnt = 0,0
            self.start_or_attach_hierarchy()
            while i < self.model_num and cnt < target:
                self.start_or_attach_hierarchy(self.hierarchy.models()[i].detached_copy())
                cnt+=1
                i += skip_num
            l.show_info('Result ensemble has {} models'.format(cnt))
            l.process_message('Renumbering models...')
            self.renum_and_count_models(hierarchy = self.new_hierarchy)
            l.process_message("Writing new ensemble to PDB file {}".format(out_name))
            self.write_hierarchy_to_pdb(out_name=out_name)
        return l
            

    def write_hierarchy_to_cif(self,
                               output_hierarchy=None,
                               symmetry=None,
                               out_name=None,
                               anisoU=False):
        """
        Write CIF file from hierarchy
        """
        from iotbx.cif import model
        cif = model.cif()
        cif_input = self.new_hierarchy.as_cif_input(crystal_symmetry=symmetry)
        cif_block = cif_input.cif_block
        cif["supercell"] = cif_block
        cif.show(out=open(out_name,'w'))
    
    def write_hierarchy_to_pdb(self,
                               output_hierarchy=None,
                               symmetry=None,
                               out_name=None,
                               anisoU=False):
        '''
        Write hierarchy to PDB file
        '''
        if symmetry == None: symmetry = self.symmetry
        if output_hierarchy == None: output_hierarchy = self.new_hierarchy
        if out_name == None: out_name = datetime_string()+'_output.pdb'
        output_hierarchy.write_pdb_file(file_name=out_name,
                                        crystal_symmetry=symmetry,
                                        anisou=anisoU)                    


class Optimiser(object):
    '''
    Simplex optimiser class. 

    Calls scitbx simplex and optimizes target function for rotation and translation
    '''
    def __init__(self,
                 trans_sigma = None,
                 rot_sigma = None,
                 trans_only = False,
                 rot_only = False,
                 pdb_object = None,
                 center_of_mass = None,
                 l = None):
        '''
        Initialize class and start simplex method
        '''
        from scitbx.simplex import simplex_opt
        # Some parameters needed in this class
        self.center_of_mass = center_of_mass
        self.pdb_object = pdb_object
        self.l = l
        self.trans_only = trans_only
        self.rot_only = rot_only
        l.process_message('Simplex minimizing rms difference B-factor profile input and ensemble')
        # Number of parameters
        n = 2
        # Create start matrix for simplex
        start_simplex = []
        if self.trans_only:
            for ii in range(n+1):
                start_simplex.append(flex.double([trans_sigma+ii*0.1,0]))
        if self.rot_only:
            for ii in range(n+1):
                start_simplex.append(flex.double([0,rot_sigma+ii*0.1]))
        if not self.trans_only and not self.rot_only:
            for ii in range(n+1):
                start_simplex.append(flex.double([trans_sigma,rot_sigma])+(ii*0.1))
        # Simplex setup
        simplex = simplex_opt(dimension = n,
                              matrix = start_simplex,
                              evaluator = self,
                              tolerance = 1.e-6)
        self.result = simplex.get_solution()
        ls = list(self.result)
        self.l.show_info('Trans_sigma = {} RMSF ({} (B-factor)), rot_sigma = {} degree'.format(ls[0],self.pdb_object.RMSF_to_B_factor(rmsf = ls[0]),ls[1]))
            
    def target(self, parameters):
        '''
        Target function simplex
        has to be called target
        '''
        trans_sigma = parameters[0]
        rot_sigma = parameters[1]
        # Call target function (PDBClass)
        score = self.pdb_object.simplex_target_func(trans_sigma = trans_sigma,
                                                    rot_sigma = rot_sigma,
                                                    center_of_mass = self.center_of_mass,
                                                    trans_only = self.trans_only,
                                                    rot_only = self.rot_only,
                                                    l = self.l)
        # Penalize negative value
        if trans_sigma < 0 or rot_sigma < 0:
            return 1e6
        # Return score
        return score

class CA_Optimiser(object):
    '''
    Simplex optimiser class. 

    Calls scitbx simplex and optimizes target function for rotation and translation
    '''
    def __init__(self,
                 trans_sigma = None,
                 rot_sigma = None,
                 pdb_object = None,
                 center_of_mass = None,
                 rb_type = None,
                 l = None):
        '''
        Initialize class and start simplex method
        '''
        from scitbx.simplex import simplex_opt
        # Some parameters needed in this class
        self.center_of_mass = center_of_mass
        self.pdb_object = pdb_object
        self.rb_type = rb_type
        self.l = l
        l.process_message('Simplex minimizing rms difference B-factor profile input and ensemble')
        # Number of parameters
        n = 2
        # Create start matrix for simplex
        start_simplex = []
        if self.rb_type == 'rot':
            for ii in range(n+1):
                start_simplex.append(flex.double([0,rot_sigma+ii*0.1]))
        if self.rb_type == 'trans':
            for ii in range(n+1):
                start_simplex.append(flex.double([trans_sigma+ii*0.1,0]))
        if self.rb_type == 'mix':
            for ii in range(n+1):
                start_simplex.append(flex.double([trans_sigma,rot_sigma])+(ii*0.1))
        # Simplex setup
        simplex = simplex_opt(dimension = n,
                              matrix = start_simplex,
                              evaluator = self,
                              tolerance = 1.e-6)
        self.result = simplex.get_solution()
        ls = list(self.result)
        self.l.show_info('Trans_sigma = {} RMSF ({} (B-factor)), rot_sigma = {} degree'.format(ls[0],self.pdb_object.RMSF_to_B_factor(rmsf = ls[0]),ls[1]))
            
    def target(self, parameters):
        '''
        Target function simplex
        has to be called target
        '''
        trans_sigma = parameters[0]
        rot_sigma = parameters[1]
        # Call target function (PDBClass)
        score = self.pdb_object.simplex_target_func_ca(trans_sigma = trans_sigma,
                                                       rot_sigma = rot_sigma,
                                                       center_of_mass = self.center_of_mass,
                                                       rb_type = self.rb_type,
                                                       l = self.l)
        # Penalize negative value
        if trans_sigma < 0 or rot_sigma < 0:
            return 1e6
        # Return score
        return score
