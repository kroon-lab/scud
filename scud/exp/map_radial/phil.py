import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    map_radial
      {
        input
          .help = "Input files"
        {
          map_in = None
            .type = path
            .help = 'map file, to do radial analysis on'
        }
        params
          .help = "Control running"
        {
          n_bins = 50
            .type = int
            .help = '# of bins in radial profile'
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
