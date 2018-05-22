import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    diffuse_all
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
          data_map = None
            .type = path
            .help = 'IDFF map created from actual data'
        }
        params
          .help = "Control running"
        {
          create_models = False
            .type = bool
            .help = 'wheter to create an ensemble from input or not'
          B_factor_plots = False
            .type = bool
            .help = 'wheter to create B-factor plots or not'
          calc_diff = False
            .type = bool
            .help = 'wheter to calculate diffuse scattering of all pdb's input or not'
          calc_CC = False
            .type = bool
            .help = 'wheter to calculate CC values between all maps and possible data or not'
          size_h = 9
            .type = int
            .help = 'number of repeats in h direction'
          size_k = 8
            .type = int
            .help = 'number of repeats in h direction'
          size_l = 5
            .type = int
            .help = 'number of repeats in h direction'
          Ncpu = 5
            .type = int
            .help = 'Number of cpus to use during parallel processes'
          supercell_num = 100
            .type = int
            .help = 'Number of supercells to create for one diffuse scattering calculation.'
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
