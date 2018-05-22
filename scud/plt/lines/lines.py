import phil as lines_phil
from scud.general.log_writer import Log
from scud.plt.Plot import Plot

def run(args=None,l=None):


    working_params = lines_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.lines

#    print p.input.dat

    if p.output.log == None:
        l = Log(log_name='lines_plot.txt')
    else:
        l = Log(log_name=p.output.log)
    l.title('plt.lines module')
    l.process_message('Processing input...\n')
    
    pt = Plot()
    for i in p.input:
        pt.read_line_data(fname=i.dat)
    if p.output.plot_out == None:
        outname = 'lines.eps'
    else:
        outname = p.output.plot_out
    pt.lines_plot(plotfile_name=outname,show=p.params.show)

    return l



