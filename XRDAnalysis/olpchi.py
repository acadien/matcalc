#!/usr/bin/python
import sys
import pylab as pl
import os,glob
import pickle
from math import asin
from numpy import *
from scipy.signal import cspline1d,cspline1d_eval
from scipy.interpolate import UnivariateSpline
from scipy.optimize import leastsq
from optparse import OptionParser

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

def localMin(seq,strt,delta):
    i = strt
    candidate = False
    m=10000.0
    base=seq[strt]
    found=False
    for elem in seq[strt:]:
        if base-elem >= delta: found=True
        if elem > base: base=elem
        if elem <= m:
            if not(candidate): base=elem
            m=elem
            candidate=True
        else:
            if found:
                return i-1
            m=elem
            candidate=False
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

def windowavg(data,N):
    return convolve(hanning(N)/sum(hanning(N)),data,mode='same')

def compute(chifile,Tfil,pressure,lamda):
    #####################
    #Get lamda and T from the temperature file
    try:
        Tfil=float(Tfil)
    except ValueError:
        fndT=False
        tlist=open(Tfil,"r")
        for line in tlist:
            if line[0]=='#':
                continue
            line=line.rstrip()
            vals=line.split(None)
            if str(vals[0])==chifile:
                Tfil=int(vals[1])
                lamda=float(vals[2])
                fndT=True
                break
        tlist.close()
        if fndT==False:
            print "Error: Unable to find chi file "+chifile+" in temperature file "+Tfil
            print "Please manually enter the temperature or enter the temperature into the temperature list file "+Tfil
            exit(0)        
    else:
        if lamda==-1.0:
            parser.print_help()
            print "\nError: When temperature is given manually, must use -l arguement to set lambda.\n"
            exit()
        if Tfil!=0 and Tfil!=300 and (Tfil<500 or Tfil>3000):
            parser.print_help()
            print "\nError: Need a valid sample temperature value for analysis"
            exit()
    
    #####################
    #Read in the chi file
    f=open(chifile,'r')
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
    i=0
    for line in f:
        line=line.lstrip().rstrip()
        [xxs[i],NULL,yys[i]]=line.split(' ')
        i+=1

    #####################
    #Load in the MgO data
    mgo_a0=4.213
    (ARs,Ts,PPs)=loadMgO('mgo_eos.dat')
    As=[i*mgo_a0 for i in ARs]
    if Tfil==0.0:
        Ps=PPs[:,0]
    elif Tfil==300.0:
        As.insert(0,mgo_a0)
        Ps=PPs[:,1].tolist()
        Ps.insert(0,1.0)
    else:
        Ps=interp_T(Tfil,Ts[2:],PPs[:,2:])
    
    #Get the MgO lattice constant for the desired pressure
    mgo_s=UnivariateSpline(Ps,As)
    mgo_sps=linspace(min(Ps),max(Ps),len(Ps))
    mgo_sas=mgo_s(mgo_sps)
    coefs=cspline1d(mgo_sas)
    latconst=cspline1d_eval(coefs,[pressure],dx=mgo_sps[1]-mgo_sps[0],x0=mgo_sps[0])

    #Get the intensity factors
    mgo_ks=[3,4,8,11,12] #Valid K vectors
    mgo_mp=[8.0,6.0,12.0,24.0,8.0] #Multiplicity
    mgo_sf=[1.0,3.0,3.0,1.0,3.0]#Structure Factor
    mgo_Im=[(mgo_sf[i]*mgo_sf[i])/mgo_mp[i] for i in range(len(mgo_sf))] #Intensity multiplicand

    #Find the MgO peaks at this pressure
    mgo_2th=list()
    mgo_I=list()
    for (i,K) in enumerate(mgo_ks):
        print lamda,K,latconst
        theta=asin(lamda*sqrt(K)*0.5/latconst)
        mgo_2th.append(degrees(theta)*2)
        mgo_I.append(mgo_Im[i]*(1+cos(theta*2)*cos(theta*2))/(sin(theta)*sin(theta)*cos(theta)))
    mxI=max(mgo_I)
    mxA=sum(xxs)/len(xxs)
    mgo_I=[(i/mxI+1)*mxA*mxA for i in mgo_I]
    pl.figure()
    pl.plot(xxs,yys)
    [pl.plot([x,x],[0,y],color='red') for (x,y) in zip(mgo_2th,mgo_I)]
    [pl.text(x,5,str(z),size='x-small') for (x,y,z) in zip(mgo_2th,mgo_I,mgo_ks)]
    pl.title(chifile+", T="+str(int(Tfil))+", P~"+str(pressure))
    pl.show()
    

##################
#Main Starts Here
##################
usage="usage: %prog [options]"
parser=OptionParser(usage=usage)
parser.add_option("-i",dest="chifilename",help="Chi input file name")
parser.add_option("-t",dest="Tfil",help="Temperature, should be T=0 or T=300 or 500<=T<=3000 or a file name where the temperature can be loaded from.")
parser.add_option("-p",dest="pressure",help="Desired pressure to plot the simulated MgO peaks")
parser.add_option("-l",dest="lamda",help="Optional, if Temperature is given manually instead of from file, wavelength of beam must be entered manually as well")

(opts,args)=parser.parse_args()

if opts.chifilename is None:
    parser.print_help()
    print "Need an input chi file from -i"
    exit(0)

if opts.Tfil is None:
    parser.print_help()
    print "\nError: Need a valid sample temperature value for proper analysis."
    exit(0)

if opts.pressure is None:
    parser.print_help()
    print "\nError: Need a pressure"
    exit(0)

if opts.lamda is None:
    opts.lamda=-1.0

compute(opts.chifilename,opts.Tfil,float(opts.pressure),float(opts.lamda))


    
