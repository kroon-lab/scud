import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    hline
      {
        input
          .help = "Input files"
        {
          mtz_in = None
            .type = path
            .help = 'File name of PDB containing ensemble to be converted to supercell'
        }
        params
          .help = "Control running"
        {
          array = None
            .type = str
            .help = 'Array name, if none, show options'
          supercell_size = 1
            .type = int
            .help = 'supercell size used for mtz creation'
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
