import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    project_map
      {
        input
          .help = "Input files"
        {
          map_in = None
            .type = path
            .help = 'input map'
        }
        params
          .help = "Control running"
        {
          mask_value = -1000
            .type = int
            .help = 'Value to mask during projection'
        }
        output
          .help = "output files"
        {
          plot_out = projection.eps
            .type = path
            .help = 'projection plot file name'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
