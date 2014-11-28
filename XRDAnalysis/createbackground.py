#!/usr/bin/python
import sys
import pylab as pl
import os
from numpy import *
from optparse import OptionParser
from scipy import interpolate
from math import *
#mine
from curvetools import *

#What does this script do:
#Given a .chi file it identifies local peaks, subtracts them from the background, smoothes the result.  
#Then writes the result to a background.chi file.

def compute():
    #####################
    #Read in the chi file
    f=open(ifilename,'r')
    i=0
    for line in f:
        if i<3:
            i+=1
            continue
        line=line.lstrip().rstrip()
        if i==3:
            i+=1
            length=int(line)
            break
    xxs=zeros(length)
    yys=zeros(length)
    i=0;
    for line in f:
        line=line.lstrip().rstrip()
        [xxs[i],NULL,yys[i]]=line.split(' ')
        i+=1
    
    s2ti= sum([1 for i in xxs if i<s2t]) #index of the starting location for 2theta

    if e2t<0:
        e2ti= len(xxs)-1
    else:
        e2ti= sum([1 for i in xxs if i<e2t]) #index of the ending location for 2theta



    #####################
    #Analyze the spectrum: find the peaks
    ind=s2ti
    inp=list()

    while True:
        ind=localMax(yys,ind,delta)
        if ind==-1 or ind>e2ti:
            break
        else:
            mxslp=max([fabs(yys[j]-yys[j+1]) for j in range(ind-5,ind+5)])
            if mxslp>0.1:
                inp.append(ind)

    rind=len(yys)-e2ti
    yreverse=yys[::-1]
    while True:
        rind=localMax(yreverse,rind,delta)
        if rind==-1 or rind<s2ti:
            break
        else:
            ind=len(yys)-rind
            found=False
            for pk in inp:
                if abs(pk-ind) < 3:
                    found=True
                    break
            if found==False:
                if ind<s2ti:
                    continue
                mxslp=max([fabs(yys[j]-yys[j+1]) for j in range(ind-5,ind+5)])
                if mxslp>0.1:
                    inp.append(ind)


    inp=sort(inp)
    peaks=[[inp[i],list(),list(),list(),list(),0] for i in range(len(inp))]
        
    #####################
    #Windowed average to smooth out the plot, useful for multi-peak analysis
    yysmth=windowavg(yys,11)

    #####################
    #Find the groups of peaks from the smooth plot
    (pstarts,pends,pgroups)=findclusterbounds(yysmth,inp)
    ingroup=list()
    for i in inp:
        found=False
        for group in pgroups:
            for j in group:
                if i==j:
                    found=True
                    break
            if found:
                break
        if not(found):
            ingroup.append(False)
        else:
            ingroup.append(True)

    ####################
    #Transfer the group peaks location from the smooth plot to the original
    for g in range(len(pgroups)):
        for s in range(len(pgroups[g])):
            smthpeak=xxs[pgroups[g][s]]
            for i in range(len(inp)):
                if abs(smthpeak-xxs[inp[i]]) < 0.05:
                    pgroups[g][s]=i
                    break

    mpeaks=list()
    ####################
    #Refine all peaks using gaussians/lorentzian fits
    i=0
    while i < len(peaks):
        ind=peaks[i][0]

        if i-1<0: lowbnd=0
        else: lowbnd=end+1

        if i+1==len(peaks): upbnd=len(xxs)-1
        else: upbnd=(peaks[i+1][0]+peaks[i][0])/2

        if min(yysmth[lowbnd:upbnd]) == 0:
            peaks.pop(i)
            inp=r_[inp[0:i],inp[i+1:len(inp)]]
            ingroup.pop(i)
            continue
        [strt,end]=findpeakbounds(yysmth,ind,lowbnd,upbnd)
        [xdata,ydata,peaks[i][4]]=glinitguess(xxs[strt:end],yys[strt:end])
        if not(ingroup[i]): #Below is an expensive function, don't call it if not necessary
            [peaks[i][2],peaks[i][3]]=glfit(xdata,ydata,peaks[i][4])
            peaks[i][1]=xdata
            peaks[i][5]=localMax(peaks[i][2],0,1)
        i=i+1

    #####################
    #Further Refine peaks by seperating out convoluted peaks
    #pl.figure()
    for i in range(len(pgroups)):
        mpeaks.append(list())
        (start,end,group)=(pstarts[i],pends[i],pgroups[i])
        groupcoefs=list()
        for p in group:
            groupcoefs.append(peaks[p][4])
        [xdata,ydata,initlsq]=mglinitguess(xxs[start:end],yys[start:end],groupcoefs)
        [mgly,lsq]=mglfit(xdata,ydata,initlsq,len(groupcoefs))
        mglx=xdata
        ilsq=extractpks(lsq,len(group))

        for j in range(len(group)):
            peaks[group[j]][3]=ilsq[j]
            [peaks[group[j]][1],peaks[group[j]][2]]=[xdata,glval(xdata,ilsq[j])]
            refpk=localMax(peaks[group[j]][2],0,1)
            peaks[group[j]][5]=refpk
        mpeaks[i].append(mglx)
        mpeaks[i].append(mgly)
            
    #####################
    #Remove these peaks from the background
    
    #First grab the index values of points that need to be removed
    ixvals=list()
    for peak in peaks:
        ixvals.append([i for i,x in enumerate(xxs) if (x>=peak[1][0] and x<=peak[1][-1])])
    for i in range(len(ixvals)):
        for j in range(8):#Grow the gaps by several points on each side
            ixvals[i].insert(0,ixvals[i][0]-1)
            ixvals[i].append(ixvals[i][-1]+1)
    ixvals= flatten(ixvals)

    #Remove bad points, creating gaps in data
    gappedis=[i for i in range(len(xxs[s2ti:e2ti])) if not(i in ixvals)]
    gappedxs=[xxs[i] for i in gappedis]
    gappedys=[yys[i] for i in gappedis]
    t=interpolate.splrep(gappedxs,gappedys,k=1) #interpolates on a line
    cappedys=interpolate.splev(xxs,t)
    
    #fit a gaussian to the background
    p0=s2ti                                        #Start point of background
    p2=finddrop(cappedys,inp[-1]+5)                #End point of background
    p1=where(cappedys==max(cappedys[p0:p2]))[0][0] #Peak of background
    smooth_cappedys=windowavg(windowavg(cappedys,20),30)
    coefs=[cappedys[p1],xxs[p1],1]
    (lsq,success)=leastsq(gaussianres,coefs,
                          args=(smooth_cappedys[p0:p2],xxs[p0:p2]),maxfev=50000)

    #append 0's where background is ignored
    background=[0 for i in xxs[:p0]] + \
               [gaussianval(x,lsq) for x in xxs[p0:p2]] + \
               [0 for i in xxs[p2:]]

    #Plot and confirm with user
    print "If you continue, the background will be written to %s"%ofilename
    print "Ctrl-C to not write this file."
    pl.figure()
    pl.plot(xxs,yys)
    pl.plot(xxs,smooth_cappedys)
    pl.plot(xxs,background)
    pl.xlabel("2Theta")
    pl.ylabel("Intensity")
    pl.title(ifilename)
    pl.legend(["Raw-Chi","Capped/Smoothed","Background"])
    pl.show()

    #############
    #Write the background to the output chifile
    import datetime
    data="This background-data chi-file was generated by %s. %s\n"%(sys.argv[0],datetime.date.today())
    data+="2-Theta Angle (Degrees)\n"
    data+="Intensity\n"
    data+="\t%d\n"%len(xxs)
    for i in range(len(xxs)):
        data+="%e  %e\n"%(xxs[i],background[i])

    ofil=open(ofilename,"w")
    ofil.write(data)


##################
#Main Starts Here
##################
usage="usage: %prog [%options]"
parser=OptionParser(usage=usage)
parser.add_option("-i",dest="ifilename",help="Input .chi file name")
parser.add_option("-o",dest="ofilename",help="Output .chi file name")
parser.add_option("-d",dest="delta",type="int",default=3,help="Optional delta arguement for peak finding, default=3")
parser.add_option("-s",dest="s2t",type="float",default=3.0,help="Starting 2Theta angle to start looking for peaks")
parser.add_option("-e",dest="e2t",type="float",default=-1,help="Ending 2Theta angle to stop looking for peaks/fitting")

(opts,args)=parser.parse_args()

if opts.ifilename is None:
    parser.print_help()
    print "Need an input chi file from -i"
    exit(0)
else:
    ifilename=opts.ifilename

if opts.ofilename is None:
    parser.print_help()
    print "Need an output chi file from -o"
    exit(0)
else:
    ofilename=opts.ofilename

if opts.delta is None:
    parser.print_help()
    print "Need a delta value for finding the peaks (restricts minimum peak height, measured from base to max)"
    exit(0)
else:
    delta=opts.delta

if opts.s2t is None:
    parser.print_help()
    print "Need a starting 2Theta angle from which to start fitting."
    exit(0)
else:
    s2t=opts.s2t

if opts.e2t is None:
    parser.print_help()
    print "Need an ending 2Theta angle from which to end fitting."
    exit(0)
else:
    e2t=opts.e2t

compute()


    
