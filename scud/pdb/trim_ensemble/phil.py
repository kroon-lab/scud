import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
      trim_ensemble
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = "File name of PDB to be trimmed"
        }
        params
          .help = "parameters for calculation"
        {
          ens_size = None
            .type = int
            .help = "number of model in final ensemble"
        }
        output
          .help = "output files"
        {
          pdb_out = "trimmed_ensemble_out.pdb"
            .type = path
            .help = 'Name of output file containing an ensemble'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
