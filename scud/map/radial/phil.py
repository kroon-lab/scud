import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    radial
      {
        input
          .help = "Input files"
        {
          map_in = None
            .type = path
            .help = 'map_file 1'
          radial_in = None
            .type = path
            .help = 'previously saved radial profile (json file)'
        }
        params
          .help = "Control running"
        {
          super_size = 5
            .type = int
            .help = 'Supercell size used in calculation or dreation of map'
          leak_num = 1
            .type = int
            .help = 'number of surrounding pixels to be corrected as wel'
          num_bins = 50
            .type = int
            .help = 'Number of bins used for radial profile'
          median = False
            .type = bool
            .help = 'Median filter from box around original Bragg'

        }
        output
          .help = "output files"
        {
          map_out = brg_filtered.map
            .type = path
            .help = 'output map name'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
