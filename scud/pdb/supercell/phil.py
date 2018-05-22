import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    supercell
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'File name of PDB containing ensemble to be converted to supercell'
          size_h = 5
            .type = int
            .help = 'Size of supercell in h direction'
          size_k = 5
            .type = int
            .help = 'Size of supercell in k direction'
          size_l = 5
            .type = int
            .help = 'Size of supercell in l direction'
          supercell_num = 5
            .type = int
            .help = 'supercell number, helps when generating loads of them'
        }
        params
          .help = "Control running"
        {
          write_pdb = False
            .type = bool
            .help = 'write PDB file?'
          write_cif = False
            .type = bool
            .help = 'write CIF file?'
          high_resolution = 2.0
            .type = float
            .help = "up to what resolution are the structurefactors calculated"
          overlap_test = False
            .type = bool
            .help = "Test for overlapping atoms"
          FFT_structure = True
            .type = bool
            .help = 'FFT supercell and write out mtz file'
          Ncpu = 30
            .type = int
            .help = 'number of cpus to use'
          max_qsub_cpu = 10
            .type = int
            .help = 'number of qsub jobs to run parallel'
        }
        output
          .help = "output files"
        {
          pdb_out = supercell_out.pdb
            .type = path
            .help = 'Name of output file containing supercell'
          mtz_out = supercell.mtz
            .type = path
            .help = "Name of output file containing structure factors"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
