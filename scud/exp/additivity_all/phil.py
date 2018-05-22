import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    additivity_all
      {
        input
          .help = "Input files"
        {
          single_pdb = None
            .type = path
            .help = 'PDB file containing single model'
          ensemble_pdb = None
            .type = path
            .help = 'PDB file containing ensemble'
        }
        params
          .help = "Control running"
        {
          create_ensemble = True
            .type = bool
            .help = 'wheter to create an ensemble from input or not'
          calc_diff = True
            .type = bool
            .help = 'Calculate diffuse scattering or not'
          merge_slices = True
            .type = bool
            .help = 'Merge 2D slices or not'
          plot_all_B_factors = True
            .type = bool
            .help = 'Plot all B-factors or not'
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
