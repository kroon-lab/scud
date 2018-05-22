import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    compare
      {
        input
          .help = "Input files"
        {
          mtz_1 = None
            .type = path
            .help = 'First mtz file'
          mtz_2 = None
            .type = path
            .help = 'Second mtz file'
        }
        params
          .help = "Control running"
        {
          array_1 = None
            .type = str
            .help = 'Array column name'
          array_2 = None
            .type = str
            .help = 'Array column name'
          supercell_size_1 = 1
            .type = int
            .help = 'Size of supercell used to generate mtz, determines miller indices that are actual Bragg'
          supercell_size_2 = 1
            .type = int
            .help = 'Size of supercell used to generate mtz, determines miller indices that are actual Bragg'
        }
        output
          .help = "output files"
        {
          mtz_1_Brg = None
            .type = path
            .help = 'when not none, write out mtz_1_Bragg positions'
          mtz_2_Brg = None
            .type = path
            .help = 'when not none, write out mtz_2_Bragg positions'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
