#!/usr/bin/python

import sys,time
from numpy import *
import scipy.interpolate as interpolate

def nPotentials(x):
    nPhi = x*(x+1)/2.
    nRho = x
    nF = x
    nU = x*(x+1)/2.
    nW = x*(x+1)/2.
    return map(int, [nPhi,nRho,nF,nU,nW])

#interpolate/extrapolate xi points from curve xx,yy
def potExtrap(xi,xx,yy):
    f=interpolate.InterpolatedUnivariateSpline(xx,yy,k=3)
    yi=f(xi)
    return yi

#read in the interpolation points for the X curves
def parse_samples(samplefile):
    sampData=open(samplefile,"r").readlines()
    datax=[[]]
    datay=[[]]
    for line in sampData:
        if len(line)<=1:
            if len(datax[-1])>0:
                datax.append(list())
                datay.append(list())
            continue
        x,y,dumb=line.split()
        datax[-1].append(float(x))
        datay[-1].append(float(y))
    datax.pop(-1)
    datay.pop(-1)
    return datax,datay

def parse_pfend(pfend):
    pfData=open(pfend,"r").readlines()

    #Parse the header
    dummy, nt, n = pfData[0].split()
    if nt!='3':
        print "This parser can only handle type 3 potfit formats for ADP potentials"
        exit(0)\

    while True:
        line = pfData.pop(0)
        if line[:2] == "#C":
            elems=line[2:].strip().split()
        if line[:2] == "#G":
            n = len(line.split())-1
            nElement = (-7+sqrt(49+24*float(n)))/6
        if line[:2] == "#E":
            break
    print "nElement types :",nElement
    nPots = nPotentials(nElement)
    xPhis,yPhis,yPhiDerivs = list(),list(),list()
    xRhos,yRhos,yRhoDerivs = list(),list(),list()
    xFs,yFs,yFDerivs = list(),list(),list()
    xUs,yUs,yUDerivs = list(),list(),list()
    xWs,yWs,yWDerivs = list(),list(),list()

    #Drop the spaces and comments
    pfData = [line for line in pfData if len(line.strip())>0 and line[0]!="#"]

    #X-Data
    genXs=lambda x:arange(x[0],x[1]+0.1,(x[1]-x[0])/(x[2]-1))
    for i in range(nPots[0]):
        l=map(float,pfData.pop(0).split())
        xPhis.append(genXs(l))

    for i in range(nPots[1]):
        l=map(float,pfData.pop(0).split())
        xRhos.append(genXs(l))

    for i in range(nPots[2]):
        l=map(float,pfData.pop(0).split())
        xFs.append(genXs(l))

    for i in range(nPots[3]):
        l=map(float,pfData.pop(0).split())
        xUs.append(genXs(l))

    for i in range(nPots[4]):
        l=map(float,pfData.pop(0).split())
        xWs.append(genXs(l))

    #Y-data
    for xPhi in xPhis:
        yPhiDerivs.append(map(float,pfData[0].split()))
        yPhis.append(array(map(float,pfData[1:size(xPhi)+1])))
        pfData=pfData[1+size(xPhi):]

    for xRho in xRhos:
        yRhoDerivs.append(map(float,pfData[0].split()))
        yRhos.append(array(map(float,pfData[1:size(xRho)+1])))
        pfData=pfData[1+size(xRho):]

    for xF in xFs:
        yFDerivs.append(map(float,pfData[0].split()))
        yFs.append(array(map(float,pfData[1:size(xF)+1])))
        pfData=pfData[1+size(xF):]

    for xU in xUs:
        yUDerivs.append(map(float,pfData[0].split()))
        yUs.append(array(map(float,pfData[1:size(xU)+1])))
        pfData=pfData[1+size(xU):]

    for xW in xWs:
        yWDerivs.append(map(float,pfData[0].split()))
        yWs.append(array(map(float,pfData[1:size(xW)+1])))
        pfData=pfData[1+size(xW):]

    #Insert new points at both ends to account for derivatives
    step=1e-3
    def appendDerivs(x,y,primes,step):
        x=insert(x,0,[x[0]-step])
        y=insert(y,0,[y[0]-primes[0]*step])

        x=insert(x,-1,[x[-1]-step])
        y=insert(y,-1,[y[-1]-primes[1]*step])
        return x,y

    xPhis,yPhis = zip(*[appendDerivs(xPhi,yPhi,yPhiDeriv,step) \
                  for xPhi,yPhi,yPhiDeriv in zip(xPhis,yPhis,yPhiDerivs)])
        
    xRhos,yRhos = zip(*[appendDerivs(xRho,yRho,yRhoDeriv,step) \
                  for xRho,yRho,yRhoDeriv in zip(xRhos,yRhos,yRhoDerivs)])

    xFs,yFs = zip(*[appendDerivs(xF,yF,yFDeriv,step) \
              for xF,yF,yFDeriv in zip(xFs,yFs,yFDerivs)])

    xUs,yUs = zip(*[appendDerivs(xU,yU,yUDeriv,step) \
              for xU,yU,yUDeriv in zip(xUs,yUs,yUDerivs)])

    xWs,yWs = zip(*[appendDerivs(xW,yW,yWDeriv,step) \
              for xW,yW,yWDeriv in zip(xWs,yWs,yWDerivs)])

    return elems,[xFs,xRhos,xPhis,xUs,xWs],[yFs,yRhos,yPhis,yUs,yWs] #return in lammps order!

def usage():
    print "%s <pf_end_??.adp> <output potential> <optional:lammps potential>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<3:
    usage()
    exit(0)

#Reformat potential and write it
elems,knotXs,knotYs=parse_pfend(sys.argv[1]) #F,rho,phi,u,w
nElem = len(elems)
cutx= max([max(map(max,i)) for i in knotXs])
Npnt=1001
dr=cutx/(Npnt-1.)
drho=1./(Npnt-1.)

axx=map(lambda x:x*dr,range(Npnt))
LRhos = [potExtrap(axx,x,y) for x,y in zip(knotXs[1],knotYs[1])]
LPhis = [potExtrap(axx,x,y) for x,y in zip(knotXs[2],knotYs[2])]
LUs = [potExtrap(axx,x,y) for x,y in zip(knotXs[3],knotYs[3])]
LWs = [potExtrap(axx,x,y) for x,y in zip(knotXs[4],knotYs[4])]

arho=[i*drho for i in range(Npnt)]
LFrhos=[potExtrap(arho,x,y) for x,y in zip(knotXs[0],knotYs[0])]

#Write element specific data
potential=[ "LAMMPS ADP potential generated by Adam Cadien\n",\
                "%s\n"%time.strftime("%d %b %Y %H:%M:%S", time.localtime()),\
                "-----\n",\
                "%d %s\n"%(nElem," ".join(elems))]
for i in range(nElem):
    potential.append("%d % 6.6f %d % 6.6f % 6.6f\n"%(Npnt,drho,Npnt,dr,cutx))
    potential.append("atomNumber mass latConst latType\n")

#LAMMPS format
for LFrho in LFrhos:
    potential+="".join(map(lambda x:"%6.6e\n"%x,LFrho))  
for LRho in LRhos:
    potential+="".join(map(lambda x:"%6.6e\n"%x,LRho))
for LPhi in LPhis:
    potential+="".join(map(lambda x:"%6.6e\n"%x,LPhi*axx))
for LU in LUs:
    potential+="".join(map(lambda x:"%6.6e\n"%x,LU))
for LW in LWs:
    potential+="".join(map(lambda x:"%6.6e\n"%x,LW))
open(sys.argv[2],"w").writelines(potential)

#If the optional arguement for the potfit lammps file is included, plots a comparison of the two
if len(sys.argv)==4:
    #below xcut points are replaced by ixpt and iypt, then interpolated over xx again.
    def cutInterp(xx,yy,xcut):
        #construct cut up arrays
        backx=xx[where(xx>xcut)[0]]
        frontx=xx[where(xx<=xcut)[0]]
        n=len(frontx)
        backy=yy[n:]
        ixpt=[0,1.0,1.5]
        y0=yy[n*2/3]*2
        y1=yy[n*2/3]
        iypt=[y0,y1,yy[n]]
        xxtemp=concatenate([ixpt,backx])
        yytemp=concatenate([iypt,backy])

        #interpolate
        f=interpolate.interp1d(xxtemp,yytemp,3)
        ynew=concatenate([f(frontx)[:n],backy])
        return ynew

    #Parse LAMMPS input Data
    potdata=open(sys.argv[3],"r").readlines()
    #Write element specific data
#    potdata[3]="1 Ge\n"
#    potdata[5]="  32 72.64 0 dummy "+" ".join(potdata[5].split()[4:])+"\n"

    nrho,drho,nr,dr,rcut = map(float,potdata[4].split())
    nr=int(nr)
    nrho=int(nrho)
    rs=array([float(i)*dr for i in range(nr)])
    rhos=array([float(i)*drho for i in range(nrho)])

    F=map(float,potdata[6:6+nrho]) #F(rho) remains unaltered
    Rho=map(float,potdata[6+nrho:6+nr+nrho])
    Phi=map(float,potdata[6+nr+nrho:6+2*nr+nrho])
    U=map(float,potdata[6+2*nr+nrho:6+3*nr+nrho])
    W=map(float,potdata[6+3*nr+nrho:6+4*nr+nrho])
    FRhoPF=array(F)
    RhoPF=array(Rho)
    PhiPF=array(Phi)/(rs+1E-16)
    UPF=array(U)
    WPF=array(W)

    #Plotting
    import pylab as pl
    pl.subplot(231)
    pl.title("Phi(r)")
    pl.scatter(knotx[2],knoty[2])
    yl=pl.ylim()
    for LPhi in LPhis:
        pl.plot(axx,LPhi,label="bleh")
    pl.plot(rs,PhiPF,label="potfit")
#    pl.ylim(yl)
    pl.legend(loc=0)

    pl.subplot(232)
    pl.title("U(r)")
    pl.scatter(knotx[3],knoty[3])
    yl=pl.ylim()
    for LU in LUs:
        pl.plot(axx,LU,label="mine")
    pl.plot(rs,UPF,label="potfit")
#    pl.ylim(yl)
    pl.legend(loc=0)

    pl.subplot(233)
    pl.title("W(r)")
    pl.scatter(knotx[4],knoty[4])
    yl=pl.ylim()
    for LW in LWs:
        pl.plot(axx,LW,label="mine")
    pl.plot(rs,WPF,label="potfit")
#    pl.ylim(yl)
    pl.legend(loc=0)

    pl.subplot(234)
    pl.title("Rho(r)")
    pl.scatter(knotx[1],knoty[1])
    yl=pl.ylim()
    for LRho in LRhos:
        pl.plot(axx,LRho,label="mine")
    pl.plot(rs,RhoPF,label="potfit")
    #pl.ylim(yl)
    pl.legend(loc=0)

    pl.subplot(235)
    pl.title("F(rho)")
    for LFrho in LFrhos:
        pl.plot(arho,LFrho,label="mine")
    pl.plot(rhos,FRhoPF,label="potfit")
    pl.scatter(knotx[0],knoty[0])
    pl.legend(loc=0)

    pl.show()

print "this code is untested for LAMMPS potentials with multiple atom types, specifically the location of the 2nd/3rd/4th chemical's information (lattice/lattice type etc) and where that goes"
