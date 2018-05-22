import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    additivity
      {
        input
          .help = "Input files"
        {
          mtz_1 = None
            .type = path
            .help = 'MTZ file, to be converted to map'
          mtz_2 = None
            .type = path
            .help = 'MTZ file, to be converted to map'
        }
        params
          .help = "Control running"
        {
          array = None
            .type = str
            .help = 'Array column name'
          sf_1 = 1.0
            .type = float
            .help = 'scale factor for array 1' 
          sf_2 = 1.0
            .type = float
            .help = 'scale factor for array 2'

        }
        output
          .help = "output files"
        {
          mtz_out = mtz_out.mtz
            .type = path
            .help = 'file name for output mtz'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
