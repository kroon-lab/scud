#!/usr/bin/env cctbx.python

from __future__ import division
import iotbx
from iotbx import pdb
import numpy as np
from scitbx.array_family import flex

class RB_PDB(object):
    '''
    Class for rigid body ensemble generation
    '''
    def __init__(self,
                 fname=None,
                 selection='All',
                 symmetry=None):

        #### Global input Parameters ####

        self.fname = fname
        self.selection = selection
        self.symmetry = symmetry
        self.n_model = None
        self.rb_hierarchy = None
        self.rb_score = np.float64(100000000.)

    def _distance_between_cart_lists(self,
                                     vec1 = None,
                                     vec2 = None):
        '''
        Calculate pair wise distances between two flex.vec3 lists
        '''
        dvec = flex.sqrt((vec1-vec2).dot())
        return dvec**2

    def _random_vector(self):
        '''
        get two random numbers
        based on that create the third
        return as np.array(vector)
        '''

        #### create random vector ####

        xy = np.random.randint(-100,high=100,size=2)/100.
        z = 1 - xy[0]**2 - xy[1]**2

        return np.array([xy[0],xy[1],z])

    def _rotation_matrix(self,
                         xyz = None,
                         angle = None):

        '''
        Create rotation matrix
        calc cos and sin first
        Then fill list in cctbx format with rotation matrix elements and return
        '''


        #### Calculate cos and Sine ####

        c = np.cos(angle)
        s = np.sin(angle)

        #### Final rotation matrix ####

        return (  c + ((xyz[0]**2)*(1-c)),
                  ((xyz[0]*xyz[1])*(1-c))-(xyz[2]*s),
                  ((xyz[0]*xyz[2])*(1-c))+(xyz[1]*s),
                  ((xyz[1]*xyz[0])*(1-c))+(xyz[2]*s),
                  c + ((xyz[1]**2)*(1-c)),
                  ((xyz[1]*xyz[2])*(1-c))-(xyz[0]*s),
                  ((xyz[2]*xyz[0])*(1-c))-(xyz[1]*s),
                  ((xyz[2]*xyz[1])*(1-c))+(xyz[0]*s),
                  c + ((xyz[2]**2)*(1-c)) )

    def calculate_center_of_mass(self):
        '''
        Calculate center of Mass from a selection/hierarchy
        '''
        self.com = self.hierarchy.atoms().extract_xyz().mean()

    def create_b_factor_mask(self,
                             n_sigma = 2):
        '''
        Mask B-factors that are n times standard deviation from the median
        '''

        self.mask = (abs(self.target_b-np.median(self.target_b)) < n_sigma * np.std(self.target_b))


    def create_rb_ens(self,
                      trans_sigma = None,
                      rot_sigma = None,
                      center_of_mass = None,
                      ens_size = None,
                      rb_type = None,
                      l = None):

        #### Create Trans and Rot lists ####

        # Rotation only
        if rb_type == 'rot':
            # zero translation vectors
            self.v_list = [[0.0,0.0,0.0] for i in range(ens_size)]
        else:
            # Random translation from multivariate normal distribution
             self.create_translation_vector_list(mu = 0,
                                                 translation_sigma = trans_sigma,
                                                 ensemble_size = ens_size)

        # rotation list
        if rb_type == 'trans':
            # zero rotation matrrices:
            rotation_list = [self._rotation_matrix(xyz = [1.0,1.0,1.0],angle = 0.0) for i in range(ens_size)]
        else:
            # Random rotations, angle from normal distribution
            rotation_list = self.create_rotation_list(rotation_sigma = rot_sigma,
                                                      ensemble_size = ens_size)

        #### Prep ensemble generation ####

        # Which models to pick in ensemble generation

        # if template is single model PDB file
        if self.n_models == 1:
            model_list = np.zeros(ens_size).astype(int)

        # If n_model is lower then desirded ensemble fill
        # with models as equally as possible, then add random models
        if self.n_models < ens_size:
            c = ens_size // self.n_models
            l = np.array(range(self.n_models))
            model_list = l
            for i in range(c-1):
                model_list = np.concatenate([model_list,l])
            rest = np.random.randint(self.n_models,
                                     size = ens_size - len(model_list))
            model_list = np.concatenate([model_list,rest])

        # if n_models is equal to ens_size just make list of model numberss
        if self.n_models == ens_size:
            model_list = np.array(range(self.n_models))

        # Create new hierarchy for rb_ensemble
        self.rb_ens_hierarchy = iotbx.pdb.hierarchy.root()

        #### Rotate/Translate Models ###

        for i,(translation_vector,rotation_matrix) in enumerate(zip(self.v_list , rotation_list)):
            # Detached model copy
            new_model = self.hierarchy.models()[model_list[i]].detached_copy()

            # rotate model
            new_model.atoms().set_xyz( ( rotation_matrix * (new_model.atoms().extract_xyz() - center_of_mass ) ) + center_of_mass )

            # translate model
            new_model.atoms().set_xyz( new_model.atoms().extract_xyz() + (translation_vector) )

            # Append to new rb_ens hierarchy
            self.rb_ens_hierarchy.append_model(new_model)


    def create_rotation_list(self,
                             rotation_sigma = None,
                             ensemble_size = None):

        '''
        Create rotation matrix list.
        Generate random angles from normal distribution.
        Then generate random vector
        Calculate rotation matrix around the vector, using angle.
        '''
        
        # Check if rotation_sigma is list or not:
        if isinstance(rotation_sigma, list) == True:

            aniso = True

            # anisotropic rotation
            rot_sigma_x = np.deg2rad(rotation_sigma[0])
            rot_sigma_y = np.deg2rad(rotation_sigma[1])
            rot_sigma_z = np.deg2rad(rotation_sigma[2])

            # Tolerance setting
            rot_sigma_tolerance = np.linalg.norm([rot_sigma_x, rot_sigma_y,rot_sigma_z])
            
        elif isinstance(rotation_sigma, list) == False:

            aniso = False

            # isotropic rotation
            rot_sigma = np.deg2rad(rotation_sigma)

            # Tolerance setting
            rot_sigma_tolerance = rot_sigma

        else:

            raise Exception('Rotation sigma not proper type')

        #### Set tolerance ####

        n = 2
        tolerance = n*abs(rot_sigma_tolerance)

        #### Prepare rotation List ####

        self.angle_list = []
        rotation_matrix_list = []

        #### Create rotation matrix list ####

        # loop over number of ensembles
        for i in range(ensemble_size):

            # Repeat if not within tolerance
            sigma_check = False

            while not sigma_check:

                # Isotropic rotation
                if aniso == False:

                    # If sigma = 0 angle is 0
                    if rot_sigma <= 1e-6:
                        angle = 0.0
                        break
                    angle = np.random.normal(0.0,rot_sigma,1)[0]

                    # Check value
                    if abs(angle) < tolerance:

                        # Escape while statement
                        sigma_check = True
   
                # anisotropic rotation
                elif aniso == True:

                    # Check for 0 rotation, then angles are 0
                    angle = np.zeros(3)
                    if rot_sigma_x >=1e-6:
                        angle[0] = np.random.normal(0.0,rot_sigma_x,1)[0]
                    if rot_sigma_y >=1e-6:
                        angle[1] = np.random.normal(0.0,rot_sigma_y,1)[0]
                    if rot_sigma_z >=1e-6:
                        angle[2] = np.random.normal(0.0,rot_sigma_z,1)[0]

                    # Check values
                    if abs(angle[0]) < tolerance and abs(angle[1]) < tolerance and abs(angle[2]) < tolerance:

                        # Escape while statement
                        sigma_check = True

                else:
                    raise Exception('Aniso boolean in rotation is wrong')

                # Random vector to rotate around:
                xyz = self._random_vector()

                #### normalize and scale if aniso ####

                if aniso == True:

                    # (1 / vectorlength) * vector
                    xyz = (1. / np.linalg.norm(xyz)) * xyz

                    # Scale angles based on vectors (sum(xyz(norm)**2)*angle)
                    angle = np.sum((xyz**2 * angle))

            # create rotation matrix from rotation vector and angle
            rotation_matrix_list.append(self._rotation_matrix(xyz,angle))

            # Save angles for plotting
            self.angle_list.append(angle)

        return rotation_matrix_list

    def create_translation_vector_list(self,
                                       mu = 0,
                                       translation_sigma = None,
                                       ensemble_size = None):

        '''
        Create a distribution of (3D) translation vectors

        - Use polar coordinates to generate xyz coordinates for translation vectors
        - Check if length of the vector is =< n*sigma
        - Depending on translation_sigma type (single value or list) do isotropic or
          anisotropic translation

        '''


        tolerance_sigma = translation_sigma

        #### Init parameters ####

        # Final vector list to return
        self.v_list = []

        # Final distance list for analysis
        self.v_dist_list = []

        # set tolerance
        n = 2 # Should be input parameter
        tolerance = n * translation_sigma 

        for i in range(ensemble_size):

            # Escape from while statement
            sigma_check = False

            # Keep repicking values untill is within tolerance
            while not sigma_check:

                #### Generate uniform random distributions ####

                v1 = np.random.random()
                v2 = 2 * np.random.random()-1
                v3 = 2 * np.random.normal()-1

                #### Factor to convert to normal distribution ####

                fac = np.sqrt(-2 * np.log(v1))

                #### Create coordinates #####

                x = (np.cos(2 * np.pi * v2) * fac) * translation_sigma
                y = (np.sin(2 * np.pi * v2) * fac) * translation_sigma
                z = (np.cos(2 * np.pi * v3) * fac) * translation_sigma

                #### Calculate distance ####

                vector_length = np.linalg.norm(np.array([x,y,z]))

                #### Check for 'extreme' values ####

                if vector_length <= tolerance:

                    # Escape while statement
                    sigma_check = True

            # Save vector and vector length to lists
            self.v_list.append(np.array([x,y,z]))
            self.v_dist_list.append(vector_length)


    def ensemble_to_B_factor(self,
                             ensemble_hierarchy=None):

        '''
        Calculate RMSF (CA) from an ensemble against the average position and convert to B-Factors
        '''

        #### Select CA's from hierarchy ####

        sel_cache = ensemble_hierarchy.atom_selection_cache()
        selection_cache = sel_cache.selection('name CA and (altid " " or altid A)')
        ensemble_hierarchy = ensemble_hierarchy.select(selection_cache)

        #### Calculate RMSF from ensemble ####

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
                    d_square = self._distance_between_cart_lists(xyz,av_xyz)
                else:
                    xyz = model.chains()[c].atoms().extract_xyz()
                    d_square += self._distance_between_cart_lists(xyz,av_xyz)
            av_d_square = d_square*(1./cnt)

            #### B from RMSF ####

            self.rb_ens_B = (8./3) * (np.pi**2) * av_d_square

        return np.array(res_num),np.array(self.rb_ens_B)

    def extract_target_b(self):
        '''
        Exctract (CA) B-factors and returns np.array
        '''

        self.target_b = np.array([at.b for at in self.hierarchy.atoms()
                                  if at.name == ' CA '])

    def read_pdb(self,
                 selection = None,
                 l = None):
        '''
        - Reads a PDB file using iotbx.pdb
        - Create hierarchy based on selection
        - Remove alt confs
        - print info

        Input:
        * filename

        Returns:
        * symmetry (global parameter)
        * pdb hierarchy (with or without selection filter)

        '''
        #### read PDB file ####

        pdb_in = iotbx.pdb.input(file_name=self.fname)

        # Read symmetry
        self.symmetry=pdb_in.crystal_symmetry()

        #### Construct hierarchy ####

        raw_hierarchy = pdb_in.construct_hierarchy()
        # Define selection and create selection cache
        sel_cache = raw_hierarchy.atom_selection_cache()
        selection_cache = sel_cache.selection(selection)

        #### Construct final hierarchy ####

        self.hierarchy = raw_hierarchy.select(selection_cache)

        #### remove alt confs ####

        self.hierarchy.remove_alt_confs(always_keep_one_conformer=True)

        #### Show Info ####

        # TODO: print more info!!!
        l.show_info('Summary for {}: '.format(self.fname))
        l.show_info('# of atoms : {}'.format(self.hierarchy.atoms().size()))

    def set_b_factor(self,
                        new_b = 0.0):
        '''
        Change B-factor of atoms in self.hierarchy
        '''

        for at in self.hierarchy.atoms():
            at.b = new_b

    def set_occupancy(self,
                      new_occ = 1.0):
        '''
        Change occupancy of atoms in self.hierarchy
        '''

        for at in self.hierarchy.atoms():
            at.occ = new_occ

    def simplex_target_func_ca(self,
                               trans_sigma = None,
                               rot_sigma = None,
                               center_of_mass = None,
                               ens_size = None,
                               target_b = None,
                               mask = None,
                               rb_type = None,
                               l = None):

        '''
        Target function for simplex minimizer,
        - create ensemble and calculate rms between original
        - B-factor profile and B-factor from rms.

        !!!! Final ensemble is saved under (self.)rb_ens_hierarchy !!!!
        '''

        #### Create ensemble ####

        self.create_rb_ens(trans_sigma = trans_sigma,
                           rot_sigma = rot_sigma,
                           center_of_mass = center_of_mass,
                           ens_size = ens_size,
                           rb_type = rb_type,
                           l = l)

        # Get b-factors
        atom,b = self.ensemble_to_B_factor(ensemble_hierarchy=self.rb_ens_hierarchy)
        # Return rmsd with target b-factors
        return np.sqrt(np.mean((target_b - b[mask])**2))


class RB_Optimiser(object):
    '''
    Simplex optimiser class.

    Calls scitbx simplex and optimizes target function for rotation and translation
    '''
    def __init__(self,
                 trans_sigma = None,
                 rot_sigma = None,
                 template_pdb = None,
                 ens_size = None,
                 target_b = None,
                 mask = None,
                 rb_type = None,
                 l = None):

        '''
        RB Simplex initialization
        '''

        from scitbx.simplex import simplex_opt

        #### Gather input parameters ####

        l.process_message('Simplex minimizing rms difference B-factor profile input and ensemble')

        # Some parameters needed in this class
        self.template_pdb = template_pdb
        self.center_of_mass = self.template_pdb.com
        self.ens_size = ens_size
        self.target_b = target_b
        self.mask = mask
        self.rb_type = rb_type
        self.l = l

        #### Initialize simplex ####

        # Number of parameters
        self.n = None

        # Create start matrix for simplex
        start_simplex = None

        l.show_info('Rigid Body motion type:')

        # Translation only
        if self.rb_type == 'trans':
            self.n = 1
            self.l.show_info('\tTranslation only')
            start_simplex = np.repeat([[trans_sigma]], self.n+1, axis=0)

        # Rotation only
        elif self.rb_type == 'rot':
            self.n = 1
            self.l.show_info('\tRotation only')
            start_simplex = np.repeat([[rot_sigma]], self.n+1, axis=0)

        # Mix
        elif self.rb_type == 'mix':
            print 'mix'
            self.n = 2
            self.l.show_info('\tTranslation and Rotation')
            start_simplex = np.repeat([[trans_sigma, rot_sigma]], self.n+1, axis=0)

        else:
            raise Exception('Invalid value')

        for ii in range(self.n):
            start_simplex[ii+1, ii] = 1.1*start_simplex[ii+1, ii]

        #### Perform Simplex minimization ####

        start_simplex = map(flex.double, start_simplex)

        simplex = simplex_opt(dimension = self.n,
                              matrix = start_simplex,
                              evaluator = self,
                              tolerance = 1.e-3)

        # Extract result
        self.result = simplex.get_solution()

    def target(self, parameters):
        '''
        Target function simplex
        has to be called target
        '''

        #### Parameters to optimize ####

        if self.n == 2:
            trans_sigma = parameters[0]
            rot_sigma = parameters[1]
        if self.n == 1 and self.rb_type == 'rot':
            rot_sigma = parameters[0]
            trans_sigma = 0.0
        if self.n == 1 and self.rb_type == 'trans':
            trans_sigma = parameters[0]
            rot_sigma = 0.0

        #### Target function ####

        #### Penalize negative value for sigma's ####
        if trans_sigma < 0 or rot_sigma < 0:
            score = 1e2

        else:
            # Target function placed in RB_PDB class
            score = self.template_pdb.simplex_target_func_ca(trans_sigma = trans_sigma,
                                                             rot_sigma = rot_sigma,
                                                             center_of_mass = self.center_of_mass,
                                                             ens_size = self.ens_size,
                                                             target_b = self.target_b,
                                                             mask = self.mask,
                                                             rb_type = self.rb_type,
                                                             l = self.l)


        self.l.show_info('| trans: {:10.4f} | rot: {:10.4f} | score: {:10.4f} |'.format(trans_sigma,
                                                                                        rot_sigma,
                                                                                        score))

        # if score is better save score and hierarchy
        if self.template_pdb.rb_score > score:
            self.template_pdb.rb_score = score
            self.template_pdb.rb_hierarchy = self.template_pdb.rb_ens_hierarchy 

        # Return score
        return score

class RB_Aniso_Optimiser(object):
    '''
    Simplex optimiser class.

    Calls scitbx simplex and optimizes target function for rotation and translation

    - Aniso in rotation (3 params)
    - optimize center of mass (3 params)
    - isotropic translation (1 param)

      7 parameters to optimize
    '''
    def __init__(self,
                 trans_sigma = None,
                 rot_sigma_array = None,
                 center_of_mass = None,
                 template_pdb = None,
                 ens_size = None,
                 target_b = None,
                 mask = None,
                 rb_type = None,
                 step = None,
                 start_step_list = None,
                 l = None):

        '''
        RB Simplex initialization
        '''

        from scitbx.simplex import simplex_opt

        #### Gather input parameters ####

        l.process_message('Simplex minimizing rms difference B-factor profile input and ensemble')

        # Some parameters needed in this class
        self.template_pdb = template_pdb
        self.center_of_mass = self.template_pdb.com
        self.ens_size = ens_size
        self.target_b = target_b
        self.mask = mask
        self.rb_type = rb_type
        self.l = l

        #### Initialize simplex ####

        # Create start matrix for simplex
        start_simplex = None
        
        # Rotation setup
        if self.rb_type == 'rot':
            # Number of parameters 
            self.n = 6
            # Combine parameters
            appended_sigma = np.append(rot_sigma_array,center_of_mass,axis=0)
            # Start array
            start_simplex = np.repeat([appended_sigma] ,self.n+1,axis=0)

            #### Determin order of magnitude of step sizes ####
            
            start_step_list = start_step_list[1:]

            # Max step size (first cycle only)
            l.show_info('Step list: {}'.format([i*step for i in start_step_list]))

        # Mix setup
        elif self.rb_type == 'mix':
            # Number of parameters 
            self.n = 7
            # Combine rot and trans start sigma (Simplex accepts single list of params)
            appended_sigma = np.append(np.array([trans_sigma]),rot_sigma_array,axis=0)
            appended_sigma = np.append(appended_sigma,center_of_mass,axis=0)
            # Start array
            start_simplex = np.repeat([appended_sigma] ,self.n+1,axis=0)

            ####  Determine order of magnitude of step sizes ###

            # Max step size (first cycle only)
            l.show_info('Step list: {}'.format([i*step for i in start_step_list]))

        else:
            raise Exception('Error in aniso_rb_opt with rb_type')

        #### Start array with step sizes applied to start values ####
        for ii in range(self.n):
            start_simplex[ii+1, ii] = start_step_list[ii]*step + start_simplex[ii+1, ii]

        # Print info
        l.process_message('Minimizing {} parameters'.format(self.n))
        l.process_message('Scaling step size with {}'.format(step))

        #### Perform Simplex minimization ####

        start_simplex = map(flex.double, start_simplex)

        simplex = simplex_opt(dimension = self.n,
                              matrix = start_simplex,
                              evaluator = self,
                              tolerance = 1e-2)

        # Extract result
        self.result = simplex.get_solution()

    def target(self, parameters):
        '''
        Target function simplex
        has to be called target
        '''

        #### Parameters to optimize ####
        
        # depending on rb_type        
        if self.rb_type == 'rot':

            # Obtain rotation parameters
            rot_sigma = list(parameters)[0:3]
            center_of_mass = list(parameters)[3:6]
            trans_sigma = 0.0

            # pentalty for negative values
            if rot_sigma[0] < 0 or rot_sigma[1] < 0 or rot_sigma[2] < 0:
                #l.warning('negative rotation sigma')
                return 1e2

        elif self.rb_type == 'mix':

            # split parameters into 3 translation sigma's and 3 rot sigma's
            trans_sigma = parameters[0]
            rot_sigma = list(parameters[1:4])
            center_of_mass = list(parameters[4:7])

            # Penalty for negative values
            if trans_sigma < 0 or rot_sigma[0] < 0 or rot_sigma[1] < 0 or rot_sigma[2] < 0:
                return 1e2


        else:
            raise Exception('Wrong rb_type in aniso simplex')

        #### Target function ####

        # Target function placed in RB_PDB class
        score = self.template_pdb.simplex_target_func_ca(trans_sigma = trans_sigma,
                                                         rot_sigma = rot_sigma,
                                                         center_of_mass = center_of_mass,
                                                         ens_size = self.ens_size,
                                                         target_b = self.target_b,
                                                         mask = self.mask,
                                                         rb_type = self.rb_type,
                                                         l = self.l)
        #### Print info per cycle ####

        if self.rb_type == 'rot':
            self.l.show_info('| Rotation: {:6.2f},{:6.2f},{:6.2f}  |  COM: {:6.2f},{:6.2f},{:6.2f}  |  score: {:6.2f} |'.format(rot_sigma[0], 
                     rot_sigma[1],
                     rot_sigma[2],
                     center_of_mass[0], 
                     center_of_mass[1], 
                     center_of_mass[2], 
                     score))

        if self.rb_type == 'mix':
            self.l.show_info('| Translation: {:6.2f}  |  Rotation: {:6.2f},{:6.2f},{:6.2f}  |  COM: {:6.2f},{:6.2f},{:6.2f}  |  score: {:6.2f} |'.format(trans_sigma, 
                    rot_sigma[0], 
                    rot_sigma[1],
                    rot_sigma[2],
                    center_of_mass[0], 
                    center_of_mass[1], 
                    center_of_mass[2], 
                    score))

        # if score is better save score and hierarchy

        if self.template_pdb.rb_score > score :
            self.template_pdb.rb_score = score
            self.template_pdb.rb_hierarchy = self.template_pdb.rb_ens_hierarchy

        # Return score
        return score
