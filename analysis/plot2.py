#!/usr/bin/python

import plotRemote as pr#mine

import sys
import pylab as pl
#mine
from datatools import windowAvg

def usage():
    print "\nUsage: %s <data-file> <x-data column><s><window size> <y-data column><s><window size>"%sys.argv[0].split("/")[-1]
    print "\nA general use plotter of 2D data. \nAttempts to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"
    print "If the column number is followed by an \"s\" then a window average is applied to that data."
    print "examples:"
    print "./plot2.py datafile 0 1"
    print "./plot2.py datafile 0s 5s50 # applies a window average of size 10 to column 0 and size 50 to column 5"
    print ""

if len(sys.argv)!=4:
    usage()
    exit(0)

fname=sys.argv[1]

#X-Data
if "s" in sys.argv[2]:
    xSmooth=True
    xCol,xWAN=sys.argv[2].split("s")
    xCol=int(xCol)
    xWAN=int(xWAN) if len(xWAN)>0 else 10
else:
    xSmooth=False
    xCol=int(sys.argv[2])

#Y-Data
if "s" in sys.argv[3]:
    ySmooth=True

    yCol,yWAN=sys.argv[3].split("s")
    yCol=int(yCol)
    yWAN=int(yWAN) if len(yWAN)>0 else 10
else:
    ySmooth=False
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

#Smoothing:
if xSmooth:
    xdata=windowAvg(xdata,xWAN)
if ySmooth:
    ydata=windowAvg(ydata,yWAN)

pl.plot(xdata,ydata)

#Try to track down column labels
labels=fraw[dataNum[0]-1].split()
if len(labels)==colN:

    if xCol!=-1:
        pl.xlabel( labels[xCol] )
        
    pl.ylabel( labels[yCol] )

pr.prshow("plot2.png")

