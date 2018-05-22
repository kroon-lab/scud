import phil as er_log_phil
from scud.general.log_writer import Log
from scud.plt.Plot import Plot
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

def run(args=None,l=None):

    if l == None:
        l = Log(log_name='er_log_plot.txt')
    l.title('plt.er_log module')

    l.process_message('Processing input...\n')
    working_params = er_log_phil.phil_parse(args=args,log=l)  # use phil to process input
    p = working_params.er_log

    # Read all lines in log file
    with open(p.input.log, 'r') as f:
        lns = f.readlines()

    # Checkpoint sentences in file
    start_line = '================= Reseting structure ensemble and total Fmodel ================\n'
    fmodel_line = '============================= Saving fmodel block =============================\n'

    # Cut to relevant part
    lns = lns[lns.index(start_line):]

    # Collect [rw,rf] per step per chunck
    cnk_ls = []
    tmp_cnk = []
    for li in lns:
        if li[0] == '~':
            splt = li.split()
            if len(splt) == 18:
                tmp_cnk.append([float(splt[11]),float(splt[12])])
        if li == fmodel_line:
            cnk_ls.append(tmp_cnk)
            tmp_cnk = []
    # Append final list as well
    cnk_ls.append(tmp_cnk)

    # Plot all
    plt.style.use('bmh')
    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.95,wspace=0.1)
    ax1 = fig.add_subplot(111)
    mpl.rcParams['axes.labelsize'] = 13
    mpl.rcParams['font.size'] = 13
    mpl.rcParams['legend.fontsize'] = 13
    ax1.legend(frameon=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
    step = int(255 / len(cnk_ls))
    for i,cnk in enumerate(cnk_ls):
        cnk = zip(*cnk)
        plt.plot(cnk[0], lw=2, ls='--', c=cm.plasma(i*step))
        plt.plot(cnk[1], lw=2, c=cm.plasma(i*step))
        plt.title(p.input.log[:-4])
        plt.ylabel('R$_{free}$')
        plt.xlabel('Timestep')
    fig.savefig(p.input.log[:-4]+'_rfactors.eps',format='eps',transparent=False)

    return l



