import libtbx.phil
from scud.general.phil_methods import init_command_line_phil

def phil_parse(args=None,log=None):
    '''
    Contains default parameters, will process commandline params and 
    changes them
    '''
    # Default parameters
    master_phil = libtbx.phil.parse("""
    ens2b
      {
        input
          .help = "Input files"
        {
          pdb_in = None
            .type = path
            .help = 'File name of PDB to be rigid body modulated'
        }
        params
          .help = "parameters for calculation"
        {
          write_av_struct = False
            .type = bool
            .help = "Wheter to write the average structure used for RMSF calc to file"
        }
        output
          .help = "output files"
        {
          plot_out = "plot.eps"
            .type = path
            .help = "name of the output plot file containing B-factors"
          plot_dat = "ens2b.json"
            .type = path
            .help = "file name for data output (json)"
        }
      }
      """)

    working_phil = init_command_line_phil(master_phil=master_phil,
                                          args=args,
                                          log=log)
    return working_phil.extract()
