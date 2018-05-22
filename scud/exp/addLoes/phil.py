import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    addLoes
      {
        input
          .help = "Input files"
        {
          single_pdb = None
            .type = path
            .help = 'Single structure PDB file'
          internal_pdb = None
            .type = path
            .help = 'PDB file containing internal motion ensemble'
          translation_pdb = trans.pdb
            .type = path
            .help = 'single translated pdb file'
          int_translation_pdb = internal_trans.pdb
            .type = path
            .help = 'ensemble translated pdb file'
        }
        params
          .help = "Control running"
        {
        do_rb = False
          .type = bool
          .help = 'Perform rb action on both single and ensemble pdbs'
        trans_sigma = 0.5
          .type = float
          .help = 'Sigma used in translation of models'
        nr_models = 250
          .type = int
          .help = 'Number of models in rb output'
        do_bfactor = False
          .type = bool
          .help = 'Calculate B-factors from the ensembles'
        do_diff_calc = False
          .type = bool
          .help = 'Calculate diffuse scattering from internal, translation and int_translation'
        size_h = 9
          .type = int
          .help = 'Super parameter in h direction'
        size_k = 8
          .type = int
          .help = 'Super parameter in k direction'
        size_l = 5
          .type = int
          .help = 'Super parameter in l direction'
        do_slice = False
          .type = bool
          .help = 'Slice diffuse scattering mtz files'
        vmin = 1000
          .type = int
          .help = 'plotting minimum'
        vmax = 75000000
          .type = int
          .help = 'plotting maximum'
        vmin_internal = 10000
          .type = int
          .help = 'plotting minimum'
        vmax_internal = 25000000
          .type = int
          .help = 'plotting maximum'
        do_subtractions = True
          .type = bool
          .help = 'subtract translation from translation+internal, and subtract internal from translation+internal'
        vmin_sub = 10000
          .type = int
          .help = 'plotting minimum after subtraction'
        vmax_sub = 40000000
          .type = int
          .help = 'plotting maximumm after subtraction'
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
