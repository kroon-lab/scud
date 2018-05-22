import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    map2mtz
      {
        input
          .help = "Input files"
        {
           map_in = None
            .type = path
            .help = 'MTZ file, to be converted to map'
        }
        params
          .help = "Control running"
        {
          make_positive = True
            .type = bool
            .help = 'If true, add minimum value to the whole map'
        }
        output
          .help = "output files"
        {
          map_out = map_out.map
            .type = path
            .help = 'file name for output map'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
