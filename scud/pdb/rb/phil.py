import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
      rb
      {
        input
          .help = "Input files"
        {
          template_pdb = None
            .type = path
            .help = 'Single structure used for rb_type er'
          target_pdb = None
            .type = path
            .help = 'PDB file from which target Bfactors will be extracted'
          rb_type = "mix"
            .type = str
            .help = 'Type of rb operations/ method to be used. options: rot, trans, mix'
        }
        params
          .help = "parameters for calculation"
        {
          filter_target_b = True
            .type = bool
            .help = 'Filter target B-factors based on deviation from mean'
          filter_n_sigma = 2
            .type = float
            .help = 'filter B-factors based on n_sigma*std'
          ensemble_size = 100
            .type = int
            .help = 'Size of output ensemble'
          seed = 123456
            .type = int
            .help = 'Seed for random number generators used in this method'
        }
        output
          .help = "output files"
        {
          pdb_out = "rb_ensemble.pdb"
            .type = path
            .help = 'Name of output file containing an ensemble'
          plot_out = "rb_ensemble.eps"
            .type = path
            .help = "name of the output plot file containing B-factors"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
