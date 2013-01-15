#!/usr/bin/python

import sys
import pylab as pl

def usage():
    print "%s <data-file> <x-data column> <y-data column>"%sys.argv[0].split("/")[-1]
    print "A general use plotter, will attempt to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"

if len(sys.argv)!=4:
    usage()
    exit(0)

data=list() 
datanums=list()
for i,line in enumerate(open(sys.argv[1],"r").readlines()):
    #drop lines with only newline characters
    if len(line) < 2: 
        continue

    #drop lines that don't contain columns of numbers
    l=line.split()
    try:
        l=map(float,l)
    except ValueError:
        continue
    data.append(l)
    datanums.append(i)

#Keep only data with the most common number of columns
lengths=map(len,data)
print [lengths.count(i) for i in range(max(lengths)+1)]
