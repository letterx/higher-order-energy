#!/usr/bin/env python

import sys

outfile = open(sys.argv[1], "w")
for fname in sys.argv[2:]:
    f = open(fname, "r")
    for line in f:
        fields = line.split()
        num_ns_edges = int(fields[10])
        num_edges = int(fields[11])
        ns_weight = int(fields[12])
        weight = int(fields[13])
        labeled = float(fields[3])/float(fields[1])

        ns_fraction = float(num_ns_edges)/float(num_edges)
        ns_weight_fraction = float(ns_weight)/float(weight)

        line = "\t".join(map(str, [labeled, ns_fraction, ns_weight_fraction]))
        outfile.write(line+'\n')
