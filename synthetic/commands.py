import subprocess
import matplotlib.pyplot as plt

# setup constants
allm = ['hocr', 'fix', 'pc']
colors = ['b^', 'rs', 'go']
cw = 2
ch = 2
alln = range(20, 400, 10)
unary = 10000000
cmd_template = './random-graph -m %s -n %d -g --unary %d --cw %d --ch %d'

energy = {}
original_variables = {}
labeled = {}
ones = {}
time = {}
additional_variables = {}
total_edges = {}

for m in allm:
    energy[m] = []
    original_variables[m] = []
    labeled[m] = []
    ones[m] = []
    time[m] = []
    additional_variables[m] = []
    total_edges[m] = []


for m in allm:
    for n in alln:
        cmd = cmd_template % (m, n, unary, cw, ch)
        out = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output = out.stdout.read().split('\n')
        energy[m].append(int(output[1].split()[1]))
        original_variables[m].append(int(output[2].split()[2]))
        labeled[m].append(int(output[3].split()[1]))
        ones[m].append(int(output[4].split()[1]))
        time[m].append(float(output[5].split()[1]))
        additional_variables[m].append(int(output[6].split()[2]))
        total_edges[m].append(int(output[7].split()[2]))
        if n % 100 == 0:
            print n

data = labeled
print data
for i in range(3):
    plt.plot(alln, data[allm[i]], colors[i])
plt.legend(allm, loc = 'upper left')
plt.axis([alln[0] - 10, alln[-1] + 10, 0, 1])
plt.autoscale(enable=True, axis='y')
plt.title('Labeled')
plt.ylabel('Labeled nodes')
plt.xlabel('Grid edge length')
plt.show()
