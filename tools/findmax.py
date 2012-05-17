#! /usr/bin/env python
"""Extract the 50 largest values of a column in a multi-column datafile.

Usage: findmax.py [index] [column] [file]

index  is the same as the gnuplot index -- blocks of data separated by double
       newlines
column is the column for which to find the 10 largest values
file   is the datafile

For example: 
$ python findmax.py 1 4

Burkhard Ritter, May 2012
"""

import re
import sys

# command line arguments
if (len(sys.argv) != 4):
    print "Usage: " + sys.argv[0] + " [index] [column] [file]"
    sys.exit(1)

selectIndex = int(sys.argv[1])
column = int(sys.argv[2])
filename = sys.argv[3]

# construct regex pattern
pattern = r"  "
for i in range(1, column):
    pattern += r"(\S+)\s+"
pattern += r"(\S+)"
regex = re.compile(pattern)

# variables
f = open(filename)
ln=0
index=0
empty=0
maxs = [(0, "") for i in range(50)]

# loop through file
for l in f:
    ln += 1
    if l == "\n":
        empty += 1
        continue
    else:
        if empty >= 2:
            index += 1
        empty = 0
    if index != selectIndex:
        continue
    if l[0] == "#":
        continue
    
    r = regex.match(l)
    if not r:
        print "Warning: regex does not match at line " + str(ln)
    v = float(r.group(column))
    minv = min(maxs, key=lambda a: a[0])
    if v > minv[0]:
        i = maxs.index(minv)
        maxs[i] = (v, l.rstrip())

# sort and print
maxs.sort(key=lambda a: a[0])
maxs.reverse()
print "# max value".ljust(15), "  line"
for i in maxs:
    print repr(i[0]).ljust(15), i[1]

f.close()
