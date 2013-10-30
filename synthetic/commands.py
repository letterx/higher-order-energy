import subprocess
import matplotlib.pyplot as plt
import collections

## setup constants
# different methods
allm = ['hocr', 'fix', 'pc']
# colors for each method
colors = ['b^', 'rs', 'go']
patterns = ['^', 's', 'o', ]
cliques = [(2,2)] #, (2,3), (2,4), (3,3)]
# range of n
alln = range(3000, 160000, 1000)
unary = 10000000
cmd_template = './random-graph -m %s -n %d --unary %d --cw %d --ch %d'

def init_dict():
    return collections.defaultdict(init_dict)

energy = init_dict()
original_variables = init_dict()
labeled = init_dict()
ones = init_dict()
time = init_dict()
additional_variables = init_dict()
total_edges = init_dict()

for m in allm:
    for c in cliques:
        energy[m][c] = []
        original_variables[m][c] = []
        labeled[m][c] = []
        ones[m][c] = []
        time[m][c] = []
        additional_variables[m][c] = []
        total_edges[m][c] = []

def initfile():
    filename = 'temp'
    write(filename)
    read(filename)

def initrun():
    for m in allm:
        for n in alln:
            for c in cliques:
                w, h = c
                cmd = cmd_template % (m, n, unary, w, h)
                out = subprocess.Popen(cmd.split(), stdout=output)
                output = out.stdout.read().split('\n')
                energy[m][c].append(int(output[1].split()[1]))
                original_variables[m][c].append(int(output[2].split()[2]))
                labeled[m][c].append(int(output[3].split()[1]))
                ones[m][c].append(int(output[4].split()[1]))
                time[m][c].append(float(output[5].split()[1]))
                additional_variables[m][c].append(int(output[6].split()[2]))
                total_edges[m][c].append(int(output[7].split()[2]))
                if n % 100 == 0:
                    print n

def write(outfile):
    scan = open(outfile, 'w')
    for m in allm:
        for n in alln:
            for c in cliques:
                w, h = c
                cmd = cmd_template % (m, n, unary, w, h)
                subprocess.Popen(cmd.split(), stdout=scan)
                if n % 100 == 0:
                    print n


def read(infile):
    scan = open(infile, 'r')
    for m in allm:
        for n in alln:
            for c in cliques:
                w, h = c
                scan.readline()
                energy[m][c].append(int(scan.readline()[1].split()[1]))
                original_variables[m][c].append(int(scan.readline()[2].split()[2]))
                labeled[m][c].append(int(scan.readline()[3].split()[1]))
                ones[m][c].append(int(scan.readline()[4].split()[1]))
                time[m][c].append(float(scan.readline()[5].split()[1]))
                additional_variables[m][c].append(int(scan.readline()[6].split()[2]))
                total_edges[m][c].append(int(scan.readline()[7].split()[2]))

def clique_size():
    data = additional_variables
    print data
    # setup plot
    ax = plt.figure().add_axes((.1, .2, .8, .7))
    xticks = [0, .33, .66, 1]
    ax.set_xticks(xticks)
    ax.set_xticklabels(cliques)
    ax.set_xlim([-.1, 1.1])
    '''
    plt.axis([alln[0] - 10, alln[-1] + 10, 0, 1])
    plt.autoscale(enable=True, axis='y')
    '''

    # plot data
    for i in range(3):
        plt.plot(xticks, [data[allm[i]][cliques[j]][0] for j in range(len(cliques))], colors[i], linestyle='-')

    # labels
    plt.title('Additional variables versus clique size for 100000 nodes')
    plt.ylabel('Number of additional variables')
    plt.xlabel('Cliques size')
    plt.legend(allm, loc = 'upper left')
    plt.show()


def number_of_nodes():
    data = time

    # setup plot
    plt.axis([alln[0] - 10, alln[-1] + 10, 0, 1])
    plt.autoscale(enable=True, axis='y')

    # plot data
    for i in range(3):
        plt.plot(alln, data[allm[i]][cliques[0]], colors[i], linestyle='-')

    # labels
    plt.title('Time versus number of nodes for clique size 2x2')
    plt.ylabel('Time (s)')
    plt.xlabel('Number of nodes')
    plt.legend(allm, loc = 'upper left')
    plt.show()

def number_of_edges():
    data = time

    print total_edges
    # plot data
    for i in range(3):
        plt.plot(total_edges[allm[i]][cliques[0]], data[allm[i]][cliques[0]], colors[i], linestyle='-')

    # labels
    plt.title('Time versus clique size for number of nodes')
    plt.ylabel('Time (s)')
    plt.xlabel('Number of nodes')
    plt.legend(allm, loc = 'upper left')
    plt.show()

initfile()
number_of_nodes()
