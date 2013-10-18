#!/usr/bin/python

import plotRemote as pr#mine

import re
import sys
import pylab as pl
#mine
from datatools import windowAvg,wsmooth
from colors import vizSpec

#returns a parsed file, keeping only rows that can be converted to floats
#keeps rows with the most common number of columns.
def parse(fname,delim=None):
    fraw = open(fname,"r").readlines()

    data=list() 
    dataIndex=list()
    for i,line in enumerate(fraw):
        #drop lines with only newline characters
        if len(line) < 2: continue
            
        #drop comment lines
        if line[0]=="#": continue
            
        #drop lines that don't contain columns of numbers
        l=line.split(delim)
        f=list()
        count=0
        for elem in l:
            try:
                f.append(float(elem))
            except ValueError:
                #f.append(elem)
                count+=1

        if count!=0:#len(l):
            continue

        data.append(f)
        dataIndex.append(i)
    
    #Keep only data with the most common number of columns
    columnLengths=map(len,data)
    colN=[columnLengths.count(i) for i in range(max(columnLengths)+1)]
    colNi=colN.index(max(colN))
    [dataNum,parsedData]=zip(*[[dataIndex[i],d] for i,d in enumerate(data) if columnLengths[i]==colNi])

    parsedData=zip(*parsedData)

    labels=fraw[dataNum[0]].split(delim)
    return labels,parsedData

def usage():
    print "\nUsage: %s <x-data column><s><window size> <y-data column><s><window size> <filenames>"%sys.argv[0].split("/")[-1]
    print "\nA general use plotter of 2D data. \nAttempts to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"
    print "If the column number is followed by an \"s\" then a window average is applied to that data."
    print "examples:"
    print "./plot2.py 0 1 datafile1 datafile2 datafile3..."
    print "./plot2.py 0s 5s50 datafile1 datafile2... # applies a window average of size 10 to column 0 and size 50 to column 5"
    print ""

if __name__=="__main__":
    pl.figure(figsize=[18,9])
    if len(sys.argv)<4:
        usage()
        exit(0)

    #X-Data
    xWAN=0
    if "s" in sys.argv[1]:
        xSmooth=True
        xCol,xWAN=sys.argv[1].split("s")
        xCol=int(xCol)
        xWAN=int(xWAN) if len(xWAN)>0 else 10
    else:
        xSmooth=False
        xCol=int(sys.argv[1])

    #Y-Data
    yWAN=0
    if "s" in sys.argv[2]:
        ySmooth=True
        yCol,yWAN=sys.argv[2].split("s")
        yCol=int(yCol)
        yWAN=int(yWAN) if len(yWAN)>0 else 10
    else:
        ySmooth=False
        yCol=int(sys.argv[2])

    #Files
    fnames=sys.argv[3:]
    try:
        fnamenumbers=map(lambda x:float("".join(re.findall('\d+',x))),fnames)
        if len(fnames) == len(fnamenumbers):
            fnames=zip(*sorted(zip(fnames,fnamenumbers),key=lambda x:x[1]))[0]
    except ValueError:
        pass

    #Colors
    colors = None
    if len(fnames)>7:
        #colors = [float2rgb_fire(i,0,len(fnames)) for i in range(len(fnames))]
        colors=True

    labels=list()
    fdatas=list()
    for fname in fnames:
        if fname[-3:]=="csv":
            l,f = parse(fname,",")
        else:
            l,f = parse(fname)

        labels.append(l)
        fdatas.append(f)

    labels=labels[0]

    #Error check on column selection
    for i,fdata in enumerate(fdatas):
        if yCol >= len(fdata):
            print "Error: Max column number is %d, but %d requested."%(len(fdata)-1,yCol)
        if xCol >= len(fdata):
            print "Error: Max column number is %d, but %d requested."%(len(fdata)-1,xCol)
        
        #Column selection
        ydata=fdata[yCol]
        if xCol==-1:
            xdata=range(len(ydata))
        else:
            xdata=fdata[xCol]

        #Smoothing:
        if xSmooth:
            xdata=windowAvg(xdata,xWAN)
        if ySmooth:
            ydata=windowAvg(ydata,yWAN)
        if xSmooth or ySmooth:
            WAN=max(xWAN,yWAN)
            xdata=xdata[WAN/2+1:WAN/-2]
            ydata=ydata[WAN/2+1:WAN/-2]

            

        #Use column labels if available
        if i==0 and len(labels)==len(fdata):
            if xCol!=-1:
                pl.xlabel( labels[xCol] )
            pl.ylabel( labels[yCol] )

        if colors==None:
            pl.plot(xdata,ydata,lw=1.5)
        else:
            pl.plot(xdata,ydata,lw=2,c=vizSpec(float(i)/len(fnames)))
    
    pl.legend(fnames,loc=0)
    pr.prshow("plot2.png")
