import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    read_header
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
          array = None
            .type = str
            .help = 'Array name, if none, show options'
          shift_origin = False
            .type = bool
            .help = 'Shift origin of map to 0,0,0 and write new map'
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
