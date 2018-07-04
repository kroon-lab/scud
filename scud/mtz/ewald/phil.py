import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    ewald
      {
        input
          .help = "Input files"
        {
          mtz_in = None
            .type = path
            .help = 'MTZ file, to be converted to map on Ewald sphere'
        }
        params
          .help = "Control running"
        {
          array = None
            .type = str
            .help = 'Array column name'
          center = None
            .type = str
            .help = 'Define center of Ewald sphere with respect to unit cell'
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
