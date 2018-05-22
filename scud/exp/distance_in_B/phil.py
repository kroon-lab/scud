import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    distance_in_B
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'PDB to check displacements'
          ref_pdb = None
            .type = path
            .help = 'PDB to calculate distances from'
        }
        params
          .help = "Control running"
        {
          tmp = None
            .type = str
            .help = 'placeholder'
        }
        output
          .help = "output files"
        {
          pdb_out = pdb_out.pdb
            .type = path
            .help = 'file name for output pdb'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
