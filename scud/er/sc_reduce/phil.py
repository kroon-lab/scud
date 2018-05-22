import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    sc_reduce
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'File name of supercell pdb (ensemble or single model)' 
          ref_pdb = None
            .type = path
            .help = 'Ref PDB to make reference model (P1) from'
          supercell_size = 2
            .type = int
            .help = 'Size of supercell (x*x*x)'
        }
        params
          .help = "Control running"
        {
          write_pdb = False
            .type = bool
            .help = 'write PDB file?'
        }
        output
          .help = "output files"
        {
          pdb_out = supercell_out.pdb
            .type = path
            .help = 'Name of output file containing supercell'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
