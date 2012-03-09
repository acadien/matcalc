#!/usr/bin/python
import sys
import pylab as pl
import os,glob
import pickle
from numpy import *
from math import asin
from scipy.signal import cspline1d,cspline1d_eval
from scipy.interpolate import UnivariateSpline
from scipy.optimize import leastsq
from optparse import OptionParser
#mine
from chitools import *

#This program gives a rough estimate of the pressure from a chi file
#without fitting any peaks or attempting to improve the accuracy
#of the estimated pressure.

#Can also load peak data from a file

Pfilename="plist.dat"
tol=0.05 #tolerance of checking if peak corresponds to MgO

def localMax(seq,strt,delta):
    i = strt
    candidate = False
    m=-10000.0
    base=seq[strt]
    found=False
    for elem in seq[strt:]:  
        if elem-base >= delta: found=True
        if elem < base: base=elem
        if elem >= m:
            if not(candidate): base=elem
            m = elem        
            candidate = True
        else:      
            if found: 
                return i-1
            m = elem        
            candidate = False 
        i+=1
    return -1    

def loadMgO(filename):
    mgof=open(filename,'r')
    for line in mgof:
        if line[0]=='#':
            continue
        break
    #First non comment line is temperatures
    line=line.rstrip()
    Ts=[float(i) for i in line.split(' ')]
    ARs=list() #a ratio (MgO FCC lattice constant ratio)
    alen=i=0
    PPs=zeros(0)
    for line in mgof:
        if line[0]=='#':
            continue
        elif alen==0:
            alen=line
            PPs=zeros((int(alen),len(Ts)))
            continue
        vals=line.split(' ')
        ARs.append(float(vals[0])**(1.0/3.0))
        PPs[i]=vals[1:]
        i+=1
    return (ARs,Ts,PPs)

def interp_T(T,Ts,PPs):
    xvals=zeros(1)
    xvals[0]=T
    Ps=list()
    for i in PPs:
        cj=cspline1d(i)
        Ps.append(cspline1d_eval(cj,xvals,dx=500.0,x0=Ts[0]))
    return Ps

def posmax(seq,key=lambda x:x):
    return max(enumerate(seq),key=lambda k:key(k[1]))[0]

def splinefit(xdata,ydata,nx):
    xs=arange(min(xdata),max(xdata),1.0/float(nx)/len(xdata))
    cj=cspline1d(ydata)
    return (xs,cspline1d_eval(cj,xs,dx=xdata[1]-xdata[0],x0=xdata[0]))

def findpeakbounds(data,pkind,lowbnd,upbnd):
    fd=0.3
    [s,e]=[max(pkind-15,0),min(pkind+16,len(data))]
    mean=(max(data[s:e])-min(data[s:e]))*0.9+min(data[s:e])
    
#Left bound
    i=pkind
    #Step over plateau
    while (fabs(data[i]-data[i-1]) < fd or data[i]>mean) and i>0:
        i-=1
    #Step down side
    while (data[i-1]-data[i]) < -1*fd and i>0:
        i-=1
    lowbnd=max([pkind-50,lowbnd])
    start=max(i,lowbnd)

#Right bound
    i=pkind
    #Step over plateau
    while (fabs(data[i]-data[i+1]) < fd or data[i]>mean) and i<len(data)-1:
        i+=1
    #Step down side
    while (data[i+1]-data[i]) < -1*fd and i<len(data)-1:   
        i+=1
    upbnd=min([pkind+50,upbnd])
    end=min([i,upbnd])
    return [start-1,end]

#Returns the index from thelist whos value is closest to theval
def findnearest(thelist,theval): 
    deltas=[fabs(aval-theval) for aval in thelist]
    return deltas.index(min(deltas))

###################################################################
#Begin Primary Compute Function
###################################################################
def compute(chifile,peaksfile,Tfil,delta,ksq,writefile,minp,maxp,lamda):
    xxs,yys=readchi(chifile)

    #####################
    #Read in the peaks file
    found=False
    peaksx=list()
    peaksy=list()
    approx_p=0.0
    if not(peaksfile is None):
        f=open(peaksfile,'r')
        for line in f:
            if line[0]=='#':
                continue
            line=line.lstrip().rstrip()
            [fname,Ttemp,lamda,approx_p,line]=line.split(None,4)
            if fname==chifile:
                lamda=float(lamda)
                approx_p=float(approx_p)
                if Tfil is None:
                    Tfil=float(Ttemp)
                found=True
                break
        if found:
            val=line.split(None,2)
            if val[0]=='0' or val[0]=='-1': #No information on peaks
                startloc=0
            else:
                startloc=float(val[0]) #get the new peak info
                [peaksx.append(float(i)) for i in val[1].split(',')]
                peaksx.sort()
                i=0
                px=peaksx[i]
                for (j,xx) in enumerate(xxs):
                    if px<xx:
                        peaksy.append(((xx-px)*yys[j]+(px-xxs[j-1])*yys[j-1])/(xxs[j]-xxs[j-1]))
                        i+=1
                        if i==len(peaksx):
                            break
                        px=peaksx[i]
        else:
            print "\nError: Unable to find chi file "+chifile+" in peaks file "+peaksfile
            print "\nTry calling again without the peaks file arguement."
            exit(0)
    else:
        startloc=xxs[0]

    #####################
    #Load in the MgO data
    mgo_a0=4.213
    (ARs,Ts,PPs)=loadMgO('/home/acadien/Documents/ascwork/EOS/mgo_eos.dat')
    As=[i*mgo_a0 for i in ARs]
    if Tfil==0.0:
        Ps=PPs[:,0]
    elif Tfil==300.0:
        As.insert(0,mgo_a0)
        Ps=PPs[:,1].tolist()
        Ps.insert(0,1.0)
    else:
        Ps=interp_T(Tfil,Ts[2:],PPs[:,2:])

    #####################
    #Analyze the spectrum: find the peaks
    inp=list()
    
    startloc=max(startloc,5.0)
    for (ind,xx) in enumerate(xxs):
        if xx>=startloc:
            break
    strtind=ind

    while True:
        ind=localMax(yys,ind,delta)
        if ind==-1:
            break
        else:
            mxslp=max([fabs(yys[j]-yys[j+1]) for j in range(ind-5,ind+5)])
            if mxslp>0.1:
                inp.append(ind)

    rind=50
    yreverse=yys[::-1]
    while True:
        rind=localMax(yreverse,rind,delta)
        if len(xxs)-rind>=strtind:
            break
        if rind==-1:
            break
        else:
            ind=len(yys)-rind
            found=False
            for pk in inp:
                if abs(pk-ind) < 3:
                    found=True
                    break
            if found==False:
                if ind<50:
                    continue
                mxslp=max([fabs(yys[j]-yys[j+1]) for j in range(ind-5,ind+5)])
                if mxslp>0.1:
                    inp.append(ind)

    inp.sort()

    #####################
    #Plot normal peaks
    #Draw the blue peaks first to be overwritten by red peaks
    pl.figure()
    pl.plot(xxs,yys,ls='dotted')
    xx=[xxs[i] for i in inp]
    yy=[yys[i] for i in inp]
    [pl.scatter(a,b) for (a,b) in zip(xx,yy)]
    [pl.scatter(a,b) for (a,b) in zip(peaksx,peaksy)]
    pl.title(chifile+", T="+str(Tfil)+", Kmag=sqrt("+str(ksq)+")")
    
    #####################
    #Calculate the MgO Peaks at this pressure
    #Get the MgO lattice constant for the desired pressure
    mgo_s=UnivariateSpline(Ps,As)
    print Ps
    mgo_sps=linspace(min(Ps),max(Ps),len(Ps))
    mgo_sas=mgo_s(mgo_sps)
    coefs=cspline1d(mgo_sas)
    latconst=cspline1d_eval(coefs,[approx_p],dx=mgo_sps[1]-mgo_sps[0],x0=mgo_sps[0])[0]

    #Get the intensity factors
    mgo_ks=[3,4,8,9,11,12] #Valid K vectors for MgO
    mgo_mp=[8.0,6.0,12.0,6.0,24.0,8.0] #Multiplicity
    mgo_sf=[1.0,3.0,3.0,1.0,1.0,3.0]#Structure Factor
    mgo_Im=[(mgo_sf[i]*mgo_sf[i])/mgo_mp[i] for i in range(len(mgo_sf))] #Intensity multiplier

    #Find the MgO peaks at this pressure
    mgo_2th=list()
    mgo_I=list()
    for (i,K) in enumerate(mgo_ks):
        theta=asin(lamda*sqrt(float(K))*0.5/latconst)
        mgo_2th.append(degrees(theta)*2)
        mgo_I.append(mgo_Im[i]*(1+cos(theta*2)*cos(theta*2))/(sin(theta)*sin(theta)*cos(theta)))
    mxI=max(mgo_I)
    mxA=sum(xxs)/len(xxs)
    mgo_I=[(i/mxI+1)*mxA*mxA for i in mgo_I]

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

    peaksx.extend([xxs[i] for i in inp])
    peaksy.extend([yys[i] for i in inp])

    for j in ksq: #for each k value given 

        thepeak=mgo_2th[mgo_ks.index(next(p for p in mgo_ks if int(j)==p))]
        i=findnearest(peaksx,thepeak)

        #####################
        #Refine the peak using a Spline fit
        ind=findnearest(xxs,peaksx[i])
        lowbnd=ind-10
        upbnd=ind+10
        #[strt,end]=findpeakbounds(yys,findnearest(xxs,peaksx[i]),lowbnd,upbnd)
        #[xd,yd]=splinefit(xxs[strt:end],yys[strt:end],10*(end-strt))
        
        #pl.plot(xd,yd)
        pl.plot([thepeak,thepeak],[0,peaksy[i]],color='red')
        pl.text(thepeak,5,str(j))
        
        #pkind=posmax(yd)
        #(theta2,intens)=(xd[pkind],yd[pkind])
        (theta2,intens)=(xxs[ind],yys[ind])
        """if round(theta2,2)==8.84: theta2=8.73"""
        
        d=lamda/(2.0*sin(radians(theta2)/2.0))
        
        a=d*sqrt(float(j))
        if (a>min(As)-tol and a<max(As)+tol):
            #Found a valid MgO peak
            pres=cspline1d_eval(coefs,[a],dx=sas[1]-sas[0],x0=sas[0])[0]
            if pres<minp or pres>maxp:
                continue
            ks.append(j)
            theas.append(a)
            thetas.append(theta2)

            pressure.append(pres)
            pl.scatter(theta2,intens,s=30,c='red',marker='o')
            pl.text(theta2,intens,str(round(theta2,2)),position=(theta2,intens),size='x-small')

    #####################
    #Print Results

    print "Possible Matches near approximated pressure "+str(approx_p)
    print "kmag\t| theta2(deg)\t| pressure(GPa)\t| a (Angstroms"
    for i in range(len(ks)):
        print str(ks[i])+"\t| "+str(round(thetas[i],5))+"  \t| "+str(round(pressure[i],2))+"   \t| "+str(round(theas[i],5))
    avgp=sum(pressure)/len(pressure)
    print "At temperature "+str(Tfil)+"K"
    err=max( [fabs(i-avgp) for i in pressure] )
    print "Average Pressure: "+str(round(avgp,2))+"GPa, Error: "+str(round(err,2))+"GPa, StdDev: "+str(round(std(pressure),2))

    pl.show()

    ######################
    #Write results to file
    if writefile is not None:
        print "Appending results to file "+writefile+"."
        file=open(writefile,'a')
        file.write(chifile+"\t"+str(Tfil)+"\t"+str(round(avgp,2))+"\t"+str(round(err,2))+"\t"+str(ks)+"\t"+str([str(round(i,2)) for i in pressure])+"\n")

##################
#Main Starts Here
##################
usage="usage: %prog [options]"
parser=OptionParser(usage=usage)
parser.add_option("-i",dest="chifilename",help="Chi input file name")
parser.add_option("-k",dest="ksq",help="k vector squared and summed (ie 111-> 3.0, 200->4.0), valid values:[3,4,8,9,11,12]")
parser.add_option("-t",dest="Tfil",help="Optional: Set the temperature manually instead of loading from peaks file. Should be T=0 or T=300 or 500<=T<=3000.")
parser.add_option("-d",dest="delta",type="float",default=1,help="Optional: delta arguement for peak finding, default=1")
parser.add_option("-p",dest="peaksfilename",help="Optional: file name of file with peaks information for melt region.")
parser.add_option("-w",dest="write",help="Optional: write the pressure results to a file, in similar format as the temperature file")
parser.add_option("-r",dest="range",help="Optional: Range of pressure value acceptable, excludes any values below the first number and above the 2nd number.  Format expected: num1,num2")
parser.add_option("-l",dest="lamda",help="Optional, if Temperature is given manually instead of from file, wavelength of beam must be entered manually as well")
#parser.add_option("--opress",dest="en_out_press",action='callback',callback=prepout, help="enables the writing of the pressure to "+Pfilename)

(opts,args)=parser.parse_args()

if opts.chifilename is None:
    parser.print_help()
    print "Need an input chi file from -i"
    exit(0)
else:
    chifile=opts.chifilename

if opts.peaksfilename is None:
    peaksfile=None
else:
    peaksfile=opts.peaksfilename

if opts.delta is None or float(opts.delta) <= 0:
    parser.print_help()
    print "\nError: Need a valid (>0) delta value for finding the peaks (restricts minimum peak height, measured from base to max)"
    exit(0)
else:
    delta=float(opts.delta)
                
if opts.ksq is None:
    print "Using default k**2 mag of 3.0 corresponding to K=[111]"
    ksq=[3.0]
else:
    ksq=opts.ksq.split(',')

writefile=opts.write

if opts.range is None:
    minp=0.0
    maxp=50.0
else:
    [minp,maxp]=[float(i) for i in opts.range.split(',')]

if opts.lamda is None:
    opts.lamda=-1.0
    parser.print_help()
    print "Must enter a valid lamda value"
    exit(0)
else:
    lamda=float(opts.lamda)

if opts.Tfil is not None:
    Tfil=float(opts.Tfil)
    if Tfil!=0 and Tfil!=300 and (Tfil<500 or Tfil>3000):
        parser.print_help()
        print "\nError: Need a valid sample temperature value for analysis"
        exit()
else:
    Tfil=None

compute(chifile,peaksfile,Tfil,delta,ksq,writefile,minp,maxp,lamda)


    
