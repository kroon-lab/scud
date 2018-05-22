import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    correlation
      {
        input
          .help = "Input files"
        {
          map_1 = None
            .type = path
            .help = 'map_file 1'
          map_2 = None
            .type = path
            .help = 'map_file 2'
        }
        params
          .help = "Control running"
        {
          use_all = False
            .type = bool
            .help = 'Use all points in both maps'
          filter_map_1 = True
            .type = bool
            .help = 'Use all points in both maps'
        }
        output
          .help = "output files"
        {
          json_out = out.json
            .type = path
            .help = 'Save plot info to this json file'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
