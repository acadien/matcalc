#!/usr/bin/python

#mine
import plotRemote as pr
#notmine
import sys
import pylab as pl

def usage():
    print "%s <sampled potential> <lammps potential>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<3:
    usage()
    exit(0)

sampledP=open(sys.argv[1],"r").readlines()
lammpsP=open(sys.argv[2],"r").readlines()

nvars=sum([1 for i in sampledP if len(i)<=1])/2
if nvars==5:
    varnames=["Phi","rho","F","u","w"]

#Load Sampled Potential
sampX=[[]]
sampY=[[]]
for i in sampledP:
    if len(i)<=1:
        if len(sampX[-1])>0:
            sampX.append(list())
            sampY.append(list())
        continue
    x,y,dumb=i.split()
    sampX[-1].append(float(x))
    sampY[-1].append(float(y))

#Load LAMMPS Potential
lammpsX=[list() for i in range(5)]
lammpsY=[list() for i in range(5)]
lammpsP=lammpsP[4:]
Nrho, drho, Nr, dr, cutoff = map(float,lammpsP[0].split()) 
print lammpsP.pop(0)
Nrho,Nr=int(Nrho),int(Nr)
lammpsP.pop(0)
print len(lammpsX)
for j in [0,1,3,4]:
    lammpsX[j]=[i*dr for i in range(Nr)] #foo(r)
lammpsX[2]=[i*drho for i in range(Nrho)] #F(rho)
print dr,drho
#F(rho)-rho(r)-phi(r)-u(r)-w(r)
lammpsY[2]=map(float,lammpsP[:Nrho])
lammpsY[1]=map(float,lammpsP[Nrho:Nrho+Nr])
lammpsY[0]=map(lambda x:float(x[1])/(lammpsX[0][x[0]]+1e-10),enumerate(lammpsP[Nrho+Nr:Nrho+2*Nr]))
lammpsY[3]=map(float,lammpsP[Nrho+2*Nr:Nrho+3*Nr])
lammpsY[4]=map(float,lammpsP[Nrho+3*Nr:Nrho+4*Nr])


print "Variable : minval - maxval"
mxy,mny=-1E10,1E10
for i,var in enumerate(varnames):
    pl.figure()
#    pl.plot(sampX[i],sampY[i],label=var)
    pl.plot(lammpsX[i],lammpsY[i])
    pl.scatter(sampX[i],sampY[i],label=var)

    if i!=2:
        pl.xlim([sampX[i][0],sampX[i][-1]])
        pl.ylim([min(sampY[i][1:])-0.5,max(sampY[i][1:])*1.1])
    else:
        pl.xlim([sampX[i][0]-0.1,sampX[i][-1]+0.1])
        pl.ylim([max(sampY[i][:])+0.5,min(sampY[i][:])*1.1])
    print "%3.3s : %4.4e - %4.4e"%(var,min(sampY[i]),max(sampY[i]))
    pl.legend()
    pr.prshow("stuff%s.png"%var)

