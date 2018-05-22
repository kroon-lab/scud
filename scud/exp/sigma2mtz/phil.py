import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    sigma2mtz
      {
        input
          .help = "Input files"
        {
          mtz_in = None
            .type = path
            .help = 'MTZ file, add sigmas'
        }
        params
          .help = "Control running"
        {
          IDFF_none = False
            .type = bool
            .help = 'Set non Bragg Intensities to 0'
          n_sigma = 100
            .type = int
            .help = 'make non-Bragg n times smaller'
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
