import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    lines
      {
        input
          .help = "Input files"
          .multiple = True
        {
          dat = None
            .type = path
            .help = 'File name of PDB containing ensemble to be converted to supercell'
        }
        params
          .help = "Control running"
        {
          title = None
            .type = str
            .help = 'Plot Title'
          dark = False
            .type = bool
            .help = 'Plot in dark color theme'
        }
        output
          .help = "output files"
        {
          plot_out = plt.eps
            .type = path
            .help = 'Name of output plot'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
