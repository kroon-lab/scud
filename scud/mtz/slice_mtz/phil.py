import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    slice_mtz
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
            .help = 'Array type, Ibrg, Itot or Idff'
          square = False
            .type = bool
            .help = 'Plot square of slice (for F)'
          vmin = 0
            .type = int
            .help = 'Minimum value in contour plot'
          vmax = 10000000
            .type = int
            .help = 'Maximum value in contour plot'
          write_slices = False
            .type = bool
            .help = 'Also output slices as a .npy file'
        }
        output
          .help = "output files"
        {
          plt_out = slice.eps
            .type = path
            .help = 'Save plot to this eps file'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
