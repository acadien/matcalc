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

        if count==len(l):
            continue

        data.append(f)
        dataIndex.append(i)
    
    #Keep only data with the most common number of columns
    columnLengths=map(len,data)

    colN=[columnLengths.count(i) for i in range(max(columnLengths)+1)]
    colNi=colN.index(max(colN))
    [dataNum,parsedData]=zip(*[[dataIndex[i],d] for i,d in enumerate(data) if columnLengths[i]==colNi])

    parsedData=zip(*parsedData)

    labels=fraw[dataNum[0]-1].split(delim)
    return labels,parsedData

def usage():
    print "\nUsage: %s <x-data column><s><window size> <y-data column><s><window size> <filenames>"%sys.argv[0].split("/")[-1]
    print "\nA general use plotter of 2D data. \nAttempts to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"
    print "If the column number is followed by an \"s\" then a window average is applied to that data."
    print "examples:"
    print "./plot2.py datafile #by default chooses columns 0 and 1"
    print "./plot2.py 0 1 datafile1 datafile2 datafile3..."
    print "./plot2.py 0s 5s50 datafile1 datafile2... # applies a window average of size 10 to column 0 and size 50 to column 5"
    print "./plot2.py 0 1 datafile1 0 2 datafile2 datafile3 #selects columns 0,1 for df1 and 0,2 for df2 and df3" 
    print ""

if __name__=="__main__":
    pl.figure(figsize=[18,9])
    if len(sys.argv)<2:
        usage()
        print "***Error***: Insuffcient Arguements"
        exit(0)

    columnIndeces=list()
    fileIndeces=list()
    columnFileCounter=list()
    for i in range(1,len(sys.argv)):
        reresult=re.search('^\d+[s]?[\d]*$',sys.argv[i])
        try:
            if reresult.group(0) == sys.argv[i]:
                columnIndeces.append(i)
                #How many files will have these columns selected
                if len(columnIndeces)%2==0:
                    columnFileCounter.append(0)
        except AttributeError: 
            fileIndeces.append(i)
            try:
                columnFileCounter[-1]+=1
            except IndexError: pass
            
    if len(columnIndeces)!=0:
        if len(columnIndeces)%2!=0 or columnIndeces[0]!=1:
            usage()
            print "***Error***: improperly formed column selection"
            exit(0)
    else:
        columnFileCounter=[len(sys.argv[1:])]

    #Go through the complex process of generating the column numbers which will be parsed out of the files selected
    xSmooths=[]
    ySmooths=[]
    xCols=[]
    yCols=[]
    xWANs=[]
    yWANs=[]
    if len(columnIndeces)==0:
        for i in fileIndeces:
            xCols.append(0)
            yCols.append(1)
            xSmooths.append(False)
            ySmooths.append(False)
            xWANs.append(0)
            yWANs.append(0)
    else:
        xcs=[c for i,c in enumerate(columnIndeces) if i%2==0]
        ycs=[c for i,c in enumerate(columnIndeces) if (i+1)%2==0]
        for i in range(len(columnFileCounter)):
            xc=sys.argv[xcs[i]]
            if "s" in xc:
                xsm=True
                r=xc.split("s")
                xc=int(r[0])
                if r[1]=="":
                    xw=10
                else:
                    xw=int(r[1])
            else:
                xsm=False
                xc=int(xc)
                xw=0

            yc=sys.argv[ycs[i]]
            if "s" in yc:
                ysm=True
                r=yc.split("s")
                yc=int(r[0])
                if r[1]=="":
                    yw=10
                else:
                    yw=int(r[1])
            else:
                ysm=False
                yc=int(yc)
                yw=0

            for j in range(columnFileCounter[i]):
                xCols.append(xc)
                xSmooths.append(xsm)
                xWANs.append(xw)
                yCols.append(yc)
                ySmooths.append(ysm)
                yWANs.append(yw)

    #Grab the file name
    fnames=[sys.argv[i] for i in fileIndeces]
    """
    #Sorting might introduce undesirable behavior so skip it
    if len(columnFileCounter)==1: #if you're only selecting 1 column then sort the file names
        try:
            fnamenumbers=map(lambda x:float("".join(re.findall('\d+',x))),fnames)
            if len(fnames) == len(fnamenumbers):
                fnames=zip(*sorted(zip(fnames,fnamenumbers),key=lambda x:x[1]))[0]
        except ValueError:
            pass
    """

    #Colors
    colors = None
    if len(fnames)>7:
        #colors = [float2rgb_fire(i,0,len(fnames)) for i in range(len(fnames))]
        colors=True

    #Load up the data and the label guesses
    labels=list()
    fdatas=list()
    for fname in fnames:
        if fname[-3:]=="csv":
            l,f = parse(fname,",")
        else:
            l,f = parse(fname)

        labels.append(l)
        fdatas.append(f)
    label=labels[0]

    for i in range(sum(columnFileCounter)):
        fdata=fdatas[i]
        xCol=xCols[i]
        yCol=yCols[i]
        xSmooth=xSmooths[i]
        ySmooth=ySmooths[i]
        xWAN=xWANs[i]
        yWAN=yWANs[i]

        #Error check on column selection
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
                pl.xlabel( label[0] )
            pl.ylabel( label[1] )

        if colors==None:
            pl.plot(xdata,ydata,lw=1.5)
        else:
            pl.plot(xdata,ydata,lw=2,c=vizSpec(float(i)/len(fnames)))
    
    pl.legend(fnames,loc=0)
    pr.prshow("plot2.png")
