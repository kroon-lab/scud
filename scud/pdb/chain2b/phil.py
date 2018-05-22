import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    chain2b
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'File of coordinate file to change B-factors from'
          cif_in = None
            .type = path
            .help = 'File of coordinate file to change B-factors from'
        }
        params
          .help = "parameters for calculation"
        {
          set_B = True
            .type = bool
            .help = "Change B-factor in file based on chain ID"
          per_model = False
            .type = bool
            .help = 'If true increment B-factor by model number, not per chain'
        }
        output
          .help = "output files"
        {
          pdb_out = None
            .type = path
            .help = "name of the output pdb with B-factor changes"
          cif_out = None
            .type = path
            .help = "name of the output pdb with B-factor changes"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
