import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    cypa_helix
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'MTZ file, to be converted to map'
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
          pdb_out = mtz_out.mtz
            .type = path
            .help = 'file name for output pdb'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
