#!/usr/bin/python

#mine
import plotRemote as pr
#notmine
import re, sys, select
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as ax3d
#mine
from datatools import windowAvg,gaussSmooth,wsmooth,superSmooth
from colors import vizSpec

#Special case for an idiot.
if "OUTCAR" in sys.argv:
    print "***"
    print "You're an idiot. Stop running this program on an OUTCAR and run it on some column data instead."
    print "***"
    exit(0)

if len(sys.argv)==1:
    if select.select([sys.stdin,],[],[],0.0)[0]:
        import fileinput
        sys.argv += fileinput.input().readline().split()


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
    n=len(labels)
    try:
        map(float,labels)
        labels=[" "]*n
    except:
        pass

    return labels,parsedData

def usage():
    print "\nUsage: %s <x-data column><s><window size> <y-data column><s><window size> <filenames>"\
        %sys.argv[0].split("/")[-1]
    print "\nUsage: %s <x-data column><x><scale> <y-data column><x><scale> <filenames>"%sys.argv[0].split("/")[-1]
    print "\nA general use plotter of 2D data. \nAttempts to find data in column format and plot the desired columns."
    print "If x-data column is -1 then range(len(y)) is used"
    print "If the column number is followed by an <s> then a window average is applied to that data."
    print "If the column number is followed by an <x> then scaling by that value is applied"
    print "examples:"
    print "./plot2.py datafile "
    print "./plot2.py 0 1 datafile1 datafile2 datafile3"
    print "./plot2.py 0 1 datafile1 0 2 datafile2 datafile3"
    print "./plot2.py 0 1s25 datafile1     #windowed average of width 25 is applied"
    print "./plot2.py 0x0.5 1x2.0 datafile #scale of 0.5 on x-axis and scale of 2.0 on y-axis"
    print "switches: -3d, -stagger, -sort, -avg, -hist, -hist#bins, -scatter, -noLeg, -saveData, -gauss, -ghost -logx -logy"
    print "switches: -alpha#val, -title <title>"
    print ""

if __name__=="__main__":

    if len(sys.argv)<2:
        usage()
        exit(0)

    columnIndeces=list()
    fileIndeces=list()
    columnFileCounter=list()

    #Pre-parse for switches
    nbins=80
    alpha=1.0 #transparency
    switches={"-3d":False,"-stagger":False,"-sort":False,"-avg":False,"-hist":False,"-scatter":False, "-noLeg":False, "-saveData":False, "-gauss":False, "-ghost":False, "-h":False,"-alpha":None,"-logx":False,"-logy":False}
    for i in range(len(sys.argv)-1,-1,-1):
        if "-hist" in sys.argv[i]: #special case hist with nbins proceding
            try:
                nbins = int(sys.argv[i].lstrip("-hist"))
                switches["-hist"]=True
                sys.argv.pop(i)
            except (ValueError,IndexError):
                pass
        elif "-alpha" in sys.argv[i]: #special case alpha
            switches["-alpha"]=True
            alpha = float(sys.argv[i].lstrip("-alpha"))
            sys.argv.pop(i)
        elif sys.argv[i] in switches.keys(): 
            switches[sys.argv[i]]=True
            sys.argv.pop(i)
        else:
            break

    if switches["-h"]:
        usage()
        exit(0)

    #Parse for column selection and file selection.
    for i in range(1,len(sys.argv)):
        reresult=re.search('^[-]?\d+[sx]?[-]?[\d]*\.?[\d]*$',sys.argv[i])
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

    xSmoothEnables=[]
    ySmoothEnables=[]
    xScaleEnables=[]
    yScaleEnables=[]
    xScales=[]
    yScales=[]
    xCols=[]
    yCols=[]
    xWANs=[]
    yWANs=[]
    if len(columnIndeces)==0:
        for i in fileIndeces:
            xCols.append(0)
            yCols.append(1)
            xSmoothEnables.append(False)
            ySmoothEnables.append(False)
            xScaleEnables.append(False)
            yScaleEnables.append(False)
            xWANs.append(0)
            yWANs.append(0)
            xScales.append(1.)
            yScales.append(1.)
    else:
        xcs=[c for i,c in enumerate(columnIndeces) if i%2==0]
        ycs=[c for i,c in enumerate(columnIndeces) if (i+1)%2==0]

        for i in range(len(columnFileCounter)):
            xc=sys.argv[xcs[i]]
            xsc=False
            xsm=False
            xw=0
            xscale=1.
            if "s" in xc:
                xsm=True
                r=xc.split("s")
                xc=int(r[0])
                if r[1]=="":
                    xw=10
                else:
                    xw=int(r[1])
            elif "x" in xc:
                xsc=True
                r=xc.split("x")
                xc=int(r[0])
                if r[1]=="":
                    xscale=1.
                else:
                    xscale=float(r[1])
            else:
                xc=int(xc)

            yc=sys.argv[ycs[i]]
            ysc=False
            ysm=False
            yw=0
            yscale=1.0
            if "s" in yc:
                ysm=True
                r=yc.split("s")
                yc=int(r[0])
                if r[1]=="":
                    yw=10
                else:
                    yw=int(r[1])
            elif "x" in yc:
                ysc=True
                r=yc.split("x")
                yc=int(r[0])
                if r[1]=="":
                    yscale=1.
                else:
                    yscale=float(r[1])
            else:
                yc=int(yc)

            for j in range(columnFileCounter[i]):
                xCols.append(xc)
                xSmoothEnables.append(xsm)
                xWANs.append(xw)
                xScaleEnables.append(xsc)
                xScales.append(xscale)
                yCols.append(yc)
                ySmoothEnables.append(ysm)
                yWANs.append(yw)
                yScaleEnables.append(ysc)
                yScales.append(yscale)

    #Grab the file name
    fnames=[sys.argv[i] for i in fileIndeces]
    if switches['-sort']:
        #Sorting might introduce undesirable behavior so skip it
        if len(columnFileCounter)==1: #if you're only selecting 1 column then sort the file names
            try:
                fnamenumbers=map(lambda x:float(".".join(re.findall('\d+',x))),fnames)
                if len(fnames) == len(fnamenumbers):
                    fnames=zip(*sorted(zip(fnames,fnamenumbers),key=lambda x:x[1]))[0]
            except ValueError:
                pass

    #Initialize Average
    initAvg=True
    count=0.

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
        print fname
    label=labels[0]

    fig=pl.figure(figsize=[18,9])
    pl.grid()
    if switches['-3d']:
        ax=fig.gca(projection='3d')
    
    for i in range(sum(columnFileCounter)):
        fdata=fdatas[i]
        xCol=xCols[i]
        yCol=yCols[i]

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
        xSmoothEnable=xSmoothEnables[i]
        ySmoothEnable=ySmoothEnables[i]
        xWAN=xWANs[i]
        yWAN=yWANs[i]
        xdataSmooth=[]
        ydataSmooth=[]

        if xSmoothEnable:
            if switches["-gauss"]:
                xdataSmooth=superSmooth(xdata,ydata,xWAN/100.0)
            else:
                xdataSmooth=windowAvg(xdata,xWAN)
        if ySmoothEnable:
            if switches["-gauss"]:
                ydataSmooth=superSmooth(xdata,ydata,yWAN/100.0)
            else:
                ydataSmooth=windowAvg(ydata,yWAN)

        #Correct for window offset, average introduces extra points that need to be chopped off
        if not switches["-gauss"] and xSmoothEnable or ySmoothEnable: 
            WAN=max(xWAN,yWAN)
            xdataSmooth=xdataSmooth[WAN/2+1:WAN/-2]
            ydataSmooth=ydataSmooth[WAN/2+1:WAN/-2]
            xdata=xdata[WAN/2+1:WAN/-2]
            ydata=ydata[WAN/2+1:WAN/-2]

        #Scaling - multiply by constant
        xScaleEnable=xScaleEnables[i]
        yScaleEnable=yScaleEnables[i]
        xScale=xScales[i]
        yScale=yScales[i]
        
        if xScaleEnable:
            xdata=[x*xScale for x in xdata]
            xdataSmooth=[x*xScale for x in xdataSmooth]
        if yScaleEnable:
            ydata=[y*yScale for y in ydata]
            ydataSmooth=[y*yScale for y in ydataSmooth]
        if switches["-logx"]:
            xdata=np.log(xdata)
        if switches["-logy"]:
            ydata=np.log(ydata)
        if switches["-stagger"]:
            m=min(ydata)
            if i==0:
                dely=(max(ydata)-min(ydata))/2.
            ydata=[y-m+i*dely for y in ydata]
            ydataSmooth=[y-m+i*dely for y in ydataSmooth]
            #if colors==None:
            #    pl.plot(xdata,ydata,lw=1.5)
            #else:
            #    pl.plot(xdata,ydata,lw=2,c=vizSpec(float(i)/len(fnames)))
            pl.tick_params(labelleft='off')

        #Plotting
        #Use column labels if available
        if switches["-avg"] and initAvg:
            initAvg=False
            if switches["-hist"]:
                mn=min(ydata)
                mx=max(ydata)
                dely=(mx-mn)/nbins
                avgx=[i*dely+mn for i in range(0,nbins+1)]
                avgy=np.zeros(nbins+1)
            else:
                avgx=xdata
                avgy=np.zeros(len(ydata))

        if i==0 and switches["-avg"]:
            ybins = range(nbins)
            mn=min(ydata)
            mx=max(ydata)

        if i==0 and len(label)==len(fdata):
            if xCol!=-1:
                pl.xlabel( label[xCol] )
            pl.ylabel( label[yCol] )
            if switches["-hist"]:
                pl.xlabel( label[yCol] )
                pl.ylabel("count")

        if switches["-3d"]:
            if colors==None:
                ax.plot(xdata,ydata,zs=i,lw=1.5)
            else:
                ax.plot(xdata,ydata,zs=i,lw=2,c=vizSpec(float(i)/len(fnames)))

        elif switches["-hist"]:
            mn=min(ydata)
            mx=max(ydata)
            dely=(mx-mn)/nbins
            ybins=[i*dely+mn for i in range(0,nbins+1)]
            yvals=np.bincount([(y-mn)/(mx-mn)*nbins for y in ydata]).tolist()
            if switches["-avg"]:
                avgy+=np.array(yvals)
                count+=1
            else:
                pl.plot(ybins,yvals)

        elif switches["-avg"]:
            if len(avgy)!=len(ydata):
                print "Not all data is the same length, unable to average lists of different lengths."
                exit(0)
            if ySmoothEnable:
                avgy+=np.array(ydataSmooth)
            else:
                avgy+=np.array(ydata)
            count+=1

        elif switches["-scatter"]:
            if xSmoothEnable:
                xdata=xdataSmooth
            if ySmoothEnable:
                ydata=ydataSmooth

            if colors==None:
                pl.scatter(xdata,ydata,lw=1.5,label=fnames[i])
            else:
                pl.scatter(xdata,ydata,lw=2,c=vizSpec(float(i)/len(fnames)),label=fnames[i])
        else:
            if colors==None:
                cp, = pl.plot([],[])
                cc = cp.get_color()
            else:
                cc=vizSpec(float(i)/len(fnames))

            if xSmoothEnable or ySmoothEnable:
                if switches["-ghost"]:
                    if colors==None:
                        cp, = pl.plot(xdata,ydata,alpha=0.4,zorder=1)
                        cc  = cp.get_color()
                    else:
                        pass
                        pl.plot(xdata,ydata,alpha=0.4,zorder=1,c=cc)

            if xSmoothEnable and ySmoothEnable:
                pl.plot(xdataSmooth,ydataSmooth,c=cc,lw=2,label=fnames[i],alpha=alpha)

            elif ySmoothEnable:
                pl.plot(xdata,ydataSmooth,c=cc,lw=2,label=fnames[i],alpha=alpha)

            elif xSmoothEnable:
                pl.plot(xdataSmooth,ydata,c=cc,lw=2,label=fnames[i],alpha=alpha)

            else:
                pl.plot(xdata,ydata,lw=1.5,c=cc,label=fnames[i],alpha=alpha)            

    if switches["-avg"]:
        avgy=[i/count for i in avgy]
        pl.plot(avgx,avgy)
        if switches["-saveData"]:
            data=label[xCol] + " " + label[yCol] + "\n"
            for x,y in zip(avgx,avgy):
                data+=str(x)+" "+str(y)+"\n"
            open("plot2.data","w").write(data)
    
    if not switches["-noLeg"]:
        pl.legend(loc=0)
    
    pr.prshow("plot2.png")
