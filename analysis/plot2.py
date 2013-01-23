#!/usr/bin/python

import plotRemote #mine

import sys
import pylab as pl

def usage():
    print "\nUsage: %s <data-file> <x-data column> <y-data column>"%sys.argv[0].split("/")[-1]
    print "\nA general use plotter of 2D data. \nAttempts to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"
    print ""

if len(sys.argv)!=4:
    usage()
    exit(0)

fname=sys.argv[1]
xCol=int(sys.argv[2])
yCol=int(sys.argv[3])

fraw = open(fname,"r").readlines()

data=list() 
dataNum=list()
for i,line in enumerate(fraw):
    #drop lines with only newline characters
    if len(line) < 2: 
        continue

    #drop comment lines
    if line[0]=="#":
        continue

    #drop lines that don't contain columns of numbers
    l=line.split()
    try:
        l=map(float,l)
    except ValueError:
        continue

    data.append(l)
    dataNum.append(i)

#Keep only data with the most common number of columns
columns=map(len,data)
colN=[columns.count(i) for i in range(max(columns)+1)]
colN=colN.index(max(colN))
[dataNum,data]=zip(*[[dataNum[i],d] for i,d in enumerate(data) if columns[i]==colN])
data=zip(*data)

#Error check on column selection
if yCol >= len(data):
    print "Error: Max column number is %d, but %d requested."%(len(data)-1,yCol)
if xCol >= len(data):
    print "Error: Max column number is %d, but %d requested."%(len(data)-1,xCol)

#Column selection
ydata=data[yCol]
if xCol==-1:
    xdata=range(len(ydata))
else:
    xdata=data[xCol]
pl.plot(xdata,ydata)

#Try to track down column labels
labels=fraw[dataNum[0]-1].split()
if len(labels)==colN:

    if xCol!=-1:
        pl.xlabel( labels[xCol] )
        
    pl.ylabel( labels[yCol] )

pl.show()

