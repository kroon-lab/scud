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
          npy_in = None
            .type = path
            .help = 'File containing 2D numpy array, output of this program'
        }
        params
          .help = "Control running"
        {
          array = 'Idff'
            .type = str
            .help = 'Array type, Ibrg, Itot or Idff'
          vmin = None
            .type = int
            .help = 'Minimum value in contour plot'
          vmax = None
            .type = int
            .help = 'Maximum value in contour plot'
        }
        output
          .help = "output files"
        {
          mtz_out = slice.mtz
            .type = path
            .help = 'Name of output mtz containing only sliced miller indices'
          npy_out = tmp.npy
            .type = path
            .help = 'Save 2D array to this numpy file'
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
