import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    slice_map
      {
        input
          .help = "Input files"
        {
          map_in = None
            .type = path
            .help = 'File name of PDB containing ensemble to be converted to supercell'
        }
        params
          .help = "Control running"
        {
          vmin = 0
            .type = int
            .help = 'Minimum value for coloring images'
          vmax = 50
            .type = int
            .help = 'Max value for coloring images'
        }
        output
          .help = "output files"
        {
          plot_out = plot.png
            .type = path
            .help = 'Save plot to file'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
