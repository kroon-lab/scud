import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    simple_ft
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'pdbfile, should be ensemble'
        }
        params
          .help = "Control running"
        {
          FFT = True
            .type = bool
            .help = 'do FFT or not'
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
