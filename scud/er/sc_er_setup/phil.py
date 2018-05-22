import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    sc_er_setup
      {
        input
          .help = "Input files"
        {
          rb_pdb_in = None
            .type = path
            .help = 'File name of single PDB (used to generate input for rigid body)'
          single_pdb_in = None
            .type = path
            .help = 'PDB file used for perfect supercell and data'
          supercell_num = 100
            .type = int
            .help = 'supercell number, helps when generating loads of them'
        }
        params
          .help = "Control running"
        {
          rb_type = "trans"
            .type = str
            .help = 'Type of rb operations/ method to be used.'
          size_h = 2
            .type = int
            .help = 'number of unit cells stacked in h direction'
          size_k = 2
            .type = int
            .help = 'number of unit cells stacked in k direction'
          size_l = 2
            .type = int
            .help = 'number of unit cells stacked in l direction'
          Ncpu = 10
            .type = int
            .help = 'number of cpu's used for supercell calculations'
          make_rb = True
            .type = bool
            .help = 'Run rigid body module'
          prep_single = True
            .type = bool
            .help = 'Prepare single_pdb'
          single_sc = True
            .type = bool
            .help = 'Calculate scattering and supercell from single_pdb'
          single_P1 = True
            .type = bool
            .help = 'Create mtz and pdb from single_pdb in P1'
          rb_sc = True
            .type = bool
            .help = 'Calculate scattering and supercell from rb_out_pdb'
          rb_P1 = True
            .type = bool
            .help = 'Create mtz and pdb from rb_out_pdb in P1'
          add_sigma = True
            .type = bool
            .help = 'Create new mtz files containing sigma columns'
          b_fact_noise = True
            .type = bool
            .help = 'Create noise between 0.0 and 0.1 in B-factor columns to fool ptls fitting procedure of ensemble refinement'
          high_resolution = 2.0
            .type = float
            .help = "up to what resolution are the structurefactors calculated"
        }
        output
          .help = "output files"
        {
          rb_pdb_out = rb_pdb_out.pdb
            .type = path
            .help = 'Name of output file containing rigid body ensemble'
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
