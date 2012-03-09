#!/usr/bin/python
import sys
import pylab as pl
import os,glob
import pickle
from numpy import *
from scipy.signal import cspline1d,cspline1d_eval
from scipy.interpolate import UnivariateSpline
from optparse import OptionParser
from multiprocessing import Process
#mine
from curvetools import *
from chitools import *

Pfilename="plist.dat" #default
tol=0.05 #tolerance of checking if peak corresponds to MgO

def readpeaks(fname):
    #Read in the user generated peaks file
    centerxs=list()
    startxs=list()
    endxs=list()
    for line in open(fname,"r").readlines()[1:]:
        [i,s,e]=[float(a) for a in line.split()]
        centerxs.append(i)
        startxs.append(s)
        endxs.append(e)
    return centerxs,startxs,endxs

def theta2index(xxs,locs):
    indeces=list()
    for loc in locs:
        indeces.append(int(sum([1 for i in xxs if i<=loc]))-1)
    return indeces

def openemacs(*fname):
    from subprocess import call
    fname="".join(fname)
    call(["emacs",fname])

def compute():
    print ""

    #####################
    #Read in the chi file
    xxs,yys=readchi(chifilename)    
    s2ti= int(sum([1 for i in xxs if i<s2t])) #index of the starting location
    
    #####################
    #Load in the MgO data
    mgo_a0=4.213
    import os
    (ARs,Ts,PPs)=loadMgO(os.getenv("HOME")+'/Documents/ascwork/EOS/mgo_eos.dat')
    As=[i*mgo_a0 for i in ARs]
    if Tval==0.0:
        Ps=PPs[:,0]
    elif Tval==300.0:
        As.insert(0,mgo_a0)
        Ps=PPs[:,1].tolist()
        Ps.insert(0,1.0)
    else:
        Ps=interp_T(Tval,Ts[2:],PPs[:,2:])

    #####################
    #Analyze the spectrum: find the peaks
    ind=s2ti
    inpleft=list()
    inpright=list()
    while True:
        ind=localMax(yys,ind,delta)
        if ind==-1:
            break
        else:
            inpleft.append(ind)
               
    rind=s2ti
    yreverse=yys[::-1]
    while True:
        rind=localMax(yreverse,rind,delta)

        if rind==-1 or rind<s2ti:
            break
        else:
            ind=len(yys)-rind-1
            inpright.append(ind)        

    ysmth=windowavg(yys,15)
    ind=s2ti
    while True:
        ind=localMax(ysmth,ind,delta)
        if ind==-1:
            break
        else:
            inpleft.append(ind)

    inp=list()
    for i in inpright:
        possibles=[j for j in inpleft if abs(i-j)<5]
        if len(possibles)>0:
            inp.append(sum(possibles)/len(possibles))

    inp=sort(list(set(inp)))

    #Write the current peak locations to a file, to be altered by user
    towrite="peak   leftbnd   rightbnd\n"
    towrite+="".join([str(xxs[i])+" \n" for i in inp])
    fname="peaks_"+chifilename.split("_")[2].split(".")[0]+".dat"
        
    #Check if the file already exists, if so don't overwrite it!
    exists=False
    try:
        open(fname,"r")
        [centerxsold,startxsold,endxsold]=readpeaks(fname)
        centers2=theta2index(xxs,centerxsold)
        starts2= theta2index(xxs,startxsold)
        ends2=   theta2index(xxs,endxsold)
        exists=True
    except IOError as e:
        print "%s File doesn't exist creating new one."%fname
        open(fname,"w").write(towrite)
    except ValueError as e:
        pass


    #Begin interactive portion
    print "List the starting and ending points of the peaks in the file %s."%fname
    print "You can add or delete peaks as needed from here."
    print "Close the graph when done."
    p=Process(target=openemacs,args=fname)
    p.start()

    #Open up a plot so the user can select the peaks
    pl.plot(xxs,yys)
    if exists:
        for i in centers2:
            pl.scatter(xxs[i],yys[i],marker="o")
        for i in starts2:
            pl.text(xxs[i],yys[i],"S")
        for i in ends2:
            pl.text(xxs[i],yys[i],"E")
        pl.legend(["diff pattern","center","start","end"])
    else:
        for i in inp:
            pl.scatter(xxs[i],yys[i])
    pl.show()
    p.join()
    print "Thanks, got it."

    #Read in the newly user generated peaks file
    [centerxs,startxs,endxs]=readpeaks(fname)
    pcenters=theta2index(xxs,centerxs)
    pstarts= theta2index(xxs,startxs)
    pends=   theta2index(xxs,endxs)
    peaks=[[pcenters[i],list(),list(),list(),list(),0] for i in range(len(pcenters))]

    #Check if there is enough space between the two points.
    badgap=False
    for i in range(len(pcenters)):
        if pends[i]-pstarts[i]<8:
            badgap=True
            print ""
            print "Error not enough space between start and end points on peak #%d at %g."%(i+1,xxs[pcenters[i]])
            print "Gap must be at least 7 points wide otherwise fit will be inaccurate."
    if badgap:exit()


        
    """
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
        """
    mpeaks=list()
    for i in range(len(peaks)):
        print "Fitting peak at %g.\n"%xxs[ind]
        ind=peaks[i][0]
        start=pstarts[i]
        end=pends[i]
        [xdata,ydata,coefs]=glinitguess(xxs[start:end],yys[start:end])
        peaks[i][4]=coefs
        [peaks[i][2],peaks[i][3]]=glfit(xdata,ydata,coefs)
        peaks[i][1]=xdata
        peaks[i][5]=localMax(peaks[i][2],0,1) #new peak
        """
        #####################
        #Further Refine peaks by seperating out convoluted peaks

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
        """

    #####################
    #Plot normal peaks
    #Draw the blue peaks first to be overwritten by red peaks
    pl.figure()
    pl.plot(xxs,yys,ls='dotted')
    xx=[xxs[i[0]] for i in peaks]
    yy=[yys[i[0]] for i in peaks]
    [pl.scatter(xx[i],yy[i],label=str(xx[i])) for i in range(len(peaks))]
    #[pl.text(xx[i],yy[i],str(round(xx[i],2)),position=(xx[i],yy[i]),size='x-small') for i in range(len(peaks))]
    pl.title(chifilename+", T="+str(Tval)+", Kmag=sqrt("+str(ksq)+")")
        #Plot fitted peaks
    for i in range(len(peaks)):
        pl.plot(peaks[i][1],peaks[i][2])
    for i in range(len(mpeaks)):
        pl.plot(mpeaks[i][0],mpeaks[i][1],ls='dashed')

    #Plot initial guesses
    #for i in range(len(peaks)):
    #    pl.plot(peaks[i][1],glval(peaks[i][1],peaks[i][4]))
    #    print peaks[i][4]
    #for i in range(len(mpeaks)):
    #    pl.plot(mpeaks[i][0],mpeaks[i][1],ls='dashed')

    #####################
    #Find possible peak matches
    s=UnivariateSpline(As,Ps)
    sas=linspace(min(As),max(As),len(As))
    sps=s(sas)
    coefs=cspline1d(sps)
    pressure=list()
    ks=list()
    theas=list()
    thetas=list()
    i=0
    for (ind,xgf,ygf,lsq,initlsq,npk) in peaks: #for each peak found
        theta2=xgf[npk]
        intens=ygf[npk]

        i+=1
        
        d=lamda/(2.0*sin(radians(theta2)/2.0))

        for j in ksq: #for each k value given 
            a=d*sqrt(float(j))
            if (a>min(As)-tol and a<max(As)+tol):
                #######################
                #Found an MgO peak
                ks.append(j)
                theas.append(a)
                thetas.append(theta2)

                pressure.append(cspline1d_eval(coefs,[a],dx=sas[1]-sas[0],x0=sas[0]))
                pl.scatter(theta2,intens,s=30,c='red',marker='o')
                pl.text(theta2,intens-0.5,str(round(theta2,5)),position=(theta2,intens),size='x-small')
                
                #ksq.remove(j)
                #break

    #####################
    #Print Results

    print "Possible Matches:"
    print "kmag\t| theta2(deg)\t| pressure(GPa)\t| a (Angstroms"
    for i in range(len(ks)):
        print str(ks[i])+"\t| "+str(round(thetas[i],5))+"  \t| "+str(round(pressure[i][0],2))+"   \t| "+str(round(theas[i],5))
    avgp=sum(pressure)/len(pressure)
    print "At temperature "+str(Tval)+"K"
    err=max( [fabs(i-avgp) for i in pressure] )
    print "Average Pressure: "+str(round(avgp,2))+"GPa, Error: "+str(round(err,3))+"GPa, StdDev: "+str(round(std(pressure),3))
    
    pl.show()

##################
#Main Starts Here
##################
usage="usage: %prog [options]"
parser=OptionParser(usage=usage)
parser.add_option("-i",dest="filename",help="Chi input file name")
parser.add_option("-l",dest="lamda",type="float",help="Wavelength of the beam corresponding to the chifile")
parser.add_option("-d",dest="delta",type="float",default=1.0,help="Optional delta arguement for peak finding, default=3")
parser.add_option("-s",dest="s2t",type="float",default=3.0,help="Optional starting 2Theta angle to start looking for peaks")
parser.add_option("-t",dest="Tval",type="float",help="Temperature, should be T=0 or T=300 or 500<=T<=3000")
parser.add_option("-k",dest="ksq",help="k vector squared and summed (ie 111-> 3.0, 200->4.0")

(opts,args)=parser.parse_args()

if opts.filename is None:
    parser.print_help()
    print "Need an input chi file from -i"
    exit(0)
else:
    chifilename=opts.filename

if opts.Tval is None:
    parser.print_help()
    print "\nError: Need a sample temperature or temperature file."
    exit()
else:
    Tval=opts.Tval

if opts.delta is None:
    parser.print_help()
    print "Need a delta value for finding the peaks (restricts minimum peak height, measured from base to max)"
    exit(0)
else:
    delta=opts.delta

if opts.s2t is None:
    parser.print_help()
    print "Need a starting 2Theta angle from which to start looking for peaks."
    exit(0)
else:
    s2t=opts.s2t
                
if opts.ksq is None:
    print "Using default k**2 mag of 3.0 corresponding to K=[111]"
    ksq=[3.0]
else:
    ksq=opts.ksq.split(',')

if opts.lamda is None:
    parser.print_help()
    print "Need a value for lambda (the corresponding beamline wavelength"
    exit(0)
else:
    lamda=opts.lamda

if chifilename in ['All','all','ALL']:
    path="./"
    k=0
    for infile in glob.glob(os.path.join(path,'*.chi')):
        k+=1
        [Null,chifilename]=infile.split('/')
        print chifilename
        compute()
        if k>=2:
            exit()
else:
    compute()


    
