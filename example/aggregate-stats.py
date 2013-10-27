#!/usr/bin/env python

import sys

total_ediff = 0
total_labeled = 0
total_lines = 0
total_time = 0
for fname in sys.argv[1:]:
    f = open(fname, "r")
    for line in f:
        total_lines += 1
        fields = line.split()
        total_ediff += float(fields[7]) - float(fields[8])
        total_labeled += float(fields[3])/3456.0
        total_time += float(fields[5])
total_ediff /= total_lines
total_labeled /= total_lines
total_time /= total_lines

print "Average energy diff: ", total_ediff
print "Average labeled: ", total_labeled
print "Average time: ", total_time
