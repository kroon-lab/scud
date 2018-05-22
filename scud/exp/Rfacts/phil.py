import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    Rfacts
      {
        input
          .help = "Input files"
        {
          pdb_1 = None
            .type = path
            .help = 'Calc FFT from this'
          mtz_2 = None
            .type = path
            .help = 'data mtz'
          er_mtz = None
            .type = path
            .help = 'output of ensemble refinement'
        }
        params
          .help = "Control running"
        {
          array_1 = FMODEL
            .type = str
            .help = 'array from pdb_1/mtz_1 to use'
          array_2 = ITOT2
            .type = str
            .help = 'Array from data mtz to use in Rfact calc'
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
