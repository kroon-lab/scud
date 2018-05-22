import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
      rigid_body
      {
        input
          .help = "Input files"
        {
         b_ref = None
            .type = path
            .help = 'PDB file containing B-factors to which to fit the new ensemble'
          b_ext = None
            .type = float
            .help = 'Single value which will be spread over all atoms, this will be the target for fitting.'
          ensemble_size = 100
            .type = int
            .help = 'Size of output ensemble'
          ensemble_in = None
            .type = path
            .help = 'Ensemble from ensembele refinementused for rb_type er'
          pdb_in = None
            .type = path
            .help = 'Single structure used for rb_type er'
        }
        params
          .help = "parameters for calculation"
        {
          
          translation_only = False
            .type = bool
            .help = "Will generate ensemble from only translation (rotation sigma=0)"
          rotation_only = False
            .type = bool
            .help = "Will generate ensemble fomr only rotation (translation sigma=0)"
          start_rotation_sigma = 3.0
            .type = float
            .help = 'Degrees to rotate'
          rb_type = "er"
            .type = str
            .help = 'Type of rb operations/ method to be used.'
          scale_t = 1.0
            .type = float
            .help = 'scale factor for trans sigma (after fitting procedure)'
          scale_r = 1.0
            .type = float
            .help = 'scale factor for rotation sigma (after fitting procedure)'
        }
        output
          .help = "output files"
        {
          pdb_out = "ensemble_out.pdb"
            .type = path
            .help = 'Name of output file containing an ensemble'
          plot_out = "plot.eps"
            .type = path
            .help = "name of the output plot file containing B-factors"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
