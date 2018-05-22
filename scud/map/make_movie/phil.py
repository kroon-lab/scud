import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    make_movie
      {
        input
          .help = "Input files"
        {
          map_1 = None
            .type = path
            .help = 'Input (first) map for movie making'
          map_2 = None
            .type = path
            .help = 'Input (first) map for movie making'
        }
        params
          .help = "Control running"
        {
          shift_origin = False
            .type = bool
            .help = 'Shift origin of map to 0,0,0 and write new map'
          log = False
            .type = bool
            .help = Plot log of I or not
        }
        output
          .help = "output files"
        {
          movie_out = movie.mp4
            .type = path
            .help = 'Save movie to file'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
