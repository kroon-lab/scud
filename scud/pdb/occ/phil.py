import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    occ
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'File name of PDB for occupancy correction'
        }
        params
          .help = "parameters for calculation"
        {
          set_occ = False
            .type = bool
            .help = "set occupancy to specific value"
          occ_val = 1.00
            .type = float
            .help = 'value to set occupancy to'
          occ_sel = all
            .type = str
            .help = 'atom selection that to change occupancies from'
          occ_model = False
            .type = bool
            .help = 'whether or not to normalize occupncy to number of models'
          set_B = True
            .type = bool
            .help = 'Set B factors to value'
          B = 0.00
            .type = float
            .help = 'B-factor value'
        }
        output
          .help = "output files"
        {
          pdb_out = occ_corrected.pdb
            .type = path
            .help = "name of the output pdb with corrected occupancies"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
