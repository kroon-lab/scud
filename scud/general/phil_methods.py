import libtbx.phil


def init_command_line_phil(master_phil=None,args=None,log=None):
    """
    Starting class for phil command line processing
    """
# Process command line args, check if given
    if args != None:
        command_line_phil = parse_phil_commandline(master_phil=master_phil,
                                                   args=args)
        working_phil = master_phil.fetch(sources=command_line_phil)
        if log != None:
            log.write(working_phil.as_str())
    else:
        master_phil.show()
        working_phil = master_phil
        if log != None:
            log.write(working_phil.as_str())
# Return parameters
    return working_phil


def parse_phil_commandline(master_phil=None,args=None,home_scope=''):
    """
    Built an argument interpreter from a master phil containing default parameters,
    and interpret commandline arguments
    """
    # Build the interpreter
    argument_interpreter = master_phil.command_line_argument_interpreter(home_scope=home_scope)
    # Loop over command line arguments
    com_list = []
    for arg in args:
        command_line_arg = argument_interpreter.process(arg=arg)
        com_list.append(command_line_arg)
    # return interpreted commandline arguments
    return com_list
