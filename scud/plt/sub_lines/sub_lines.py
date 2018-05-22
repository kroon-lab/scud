import phil as lines_phil
from scud.general.log_writer import Log
from scud.plt.Plot import Plot

def run(args):

    l = Log(log_name='lines_plot.txt')
    l.title('plt.lines module')

    l.process_message('Processing input...\n')
    working_params = lines_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.lines
    
    pt = Plot(dark=p.params.dark)
    for i in p.input:
        pt.read_line_data(fname=i.dat)
    pt.sub_lines_plot(plotfile_name='sub_tst.eps')




