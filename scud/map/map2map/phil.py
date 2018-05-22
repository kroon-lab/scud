import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    map2map
      {
        input
          .help = "Input files"
        {
          map_1 = None
            .type = path
            .help = 'first map'
          map_2 = None
            .type = path
            .help = 'second map'
        }
        params
          .help = "Control running"
        {
          scale_maps = True
            .type = str
            .help = 'Scale map with smallest intensity to map with smalles intensity'
          ignore_exp = -1000
            .type = int
            .help = 'Value to ignore when scaling maps'
          map_1_exp = False
            .type = bool
            .help = 'Is map_1 an experimental map or note'
          map_2_exp = False
            .type = bool
            .help = 'Is map_2 an experimental map or note'
          save_all_maps = True
            .type = bool
            .help = 'Save all maps generated during program'
          make_diff_maps = True
            .type = bool
            .help = 'Create and save difference maps?'
          bragg_cc = True
            .type = bool
            .help = 'To calculate a Bragg CC value'
          radial_bins = 20
            .type = int
            .help = 'amount of bins used for making radial profile'
          size_h = 9
            .type = int
            .help = 'Number of times unit cells were stacked along the "a" direction '
          size_k = 8
            .type = int
            .help = 'Number of times unit cells were stacked along the "b" direction '
          size_l = 5
            .type = int
            .help = 'Number of times unit cells were stacked along the "c" direction '
        }
        output
          .help = "output files"
        {
          map_out = modified_map.map
            .type = path
            .help = 'Save new map'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
