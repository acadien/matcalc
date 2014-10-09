#!/usr/bin/python

import sys
from numpy import *
#mine
import plotRemote as pr

rmsdFile = sys.argv[1]
if rmsdFile[-4:]!="rmsd": exit(0)

rcut=1.5

rmsdAvg=list()
rmsdPerAtom=list()

jumpCount50=[0]*50
jumpCount100=[0]*100
jumpCount200=[0]*200
jumpCount300=[0]*300
jumpCount400=[0]*400
jumpCount500=[0]*500
jumpCount1000=[0]*1000

jumpedFlag50=[] #how many times has the atom jumped >2.5A in 50ts
jumpedFlag100=[] #how many times has the atom jumped >2.5A in 100ts
jumpedFlag200=[]
jumpedFlag300=[]
jumpedFlag400=[]
jumpedFlag500=[]
jumpedFlag1000=[]

for i,line in enumerate(open(rmsdFile,"r").readlines()):
    if i==0:
        continue
    line = map(float,line.split())
    rmsdAvg.append(line[0])
    rmsdPerAtom.append(line[1:])

    nIon=len(line[1:])
    if i==1:
        jumpedFlag50 = zeros(nIon)
        jumpedFlag100 = zeros(nIon)
        jumpedFlag200 = zeros(nIon)
        jumpedFlag300 = zeros(nIon)
        jumpedFlag400 = zeros(nIon)
        jumpedFlag500 = zeros(nIon)
        jumpedFlag1000 = zeros(nIon)

    if i>50:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-50]) if fabs(b-a)>rcut])
        jumpCount50.append(c)
        jumpedFlag50 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-50]) ])
                
    if i>100:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-100]) if fabs(b-a)>rcut])
        jumpCount100.append(c)
        jumpedFlag100 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-100]) ])

    if i>200:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-200]) if fabs(b-a)>rcut])
        jumpCount200.append(c)
        jumpedFlag200 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-200]) ])

    if i>300:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-300]) if fabs(b-a)>rcut])
        jumpCount300.append(c)
        jumpedFlag300 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-300]) ])

    if i>400:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-400]) if fabs(b-a)>rcut])
        jumpCount400.append(c)
        jumpedFlag400 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-400]) ])

    if i>500:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-500]) if fabs(b-a)>rcut])
        jumpCount500.append(c)
        jumpedFlag500 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-500]) ])

    if i>1000:
        c=sum([1 for a,b in zip(line[1:],rmsdPerAtom[i-1000]) if fabs(b-a)>rcut])
        jumpCount1000.append(c)
        jumpedFlag1000 += array([1 if fabs(b-a)>rcut else 0 for a,b, in zip(line[1:],rmsdPerAtom[i-1000]) ])


#do the fast/slow demarcation (total)
rmsdPerAtom=array(rmsdPerAtom)
slowRMSDAvg=zeros(len(rmsdAvg))
fastRMSDAvg=zeros(len(rmsdAvg))
fastAtoms=[]
fc,sc=0,0
f=open("hoppedAtomsOCFULL.dat","w")
f.write("iter j50 j100 j200 j500 j1000\n")
for i in range(nIon):
    f.write("%d %d %d %d %d %d\n"%(i,jumpedFlag50[i],jumpedFlag100[i],jumpedFlag200[i],jumpedFlag500[i],jumpedFlag1000[i]))
    if jumpedFlag50[i]>0:
        fc+=1
        fastRMSDAvg += rmsdPerAtom[:,i]
        #fastAtoms.append(rmsdPerAtom[:,i])
    else:
        sc+=1
        slowRMSDAvg += rmsdPerAtom[:,i]
f.close()
fastRMSDAvg/=fc
slowRMSDAvg/=sc

import pylab as pl
pl.figure(figsize=(20,14))

#pl.subplot(311)
pl.plot(jumpCount100,label="Jts100")    
#pl.plot(jumpCount200,label="Jts200")    
#pl.plot(jumpCount300,label="Jts300")    
pl.plot(jumpCount400,label="Jts400")    
pl.plot(jumpCount500,label="Jts500")    
pl.plot(jumpCount1000,label="Jts500")    
#pl.plot(jumpCount100,label="Jts100")    
pl.plot(rmsdAvg,lw=2,label="RMSD")
pl.legend(loc=0)
pl.ylabel("# of atoms that have Jumped >%1.1fA"%rcut)
pl.xlabel("timestep")
"""
pl.subplot(312)

pl.plot(rmsdAvg,label="AvgRMSD")
pl.plot(slowRMSDAvg,label="SlowAvgRMSD")
pl.plot(fastRMSDAvg,label="FastAvgRMSD")
for i in fastAtoms:
    pl.plot(i)
pl.legend(loc=0)

pl.subplot(313)
pl.scatter(range(nIon),jumpedFlag10)
pl.xlim([0,nIon])
"""
nJumped = sum([1 for i in jumpedFlag50 if i>0])
print "%d: %d / %d"%(50,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag100 if i>0])
print "%d: %d / %d"%(100,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag200 if i>0])
print "%d: %d / %d"%(200,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag300 if i>0])
print "%d: %d / %d"%(300,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag400 if i>0])
print "%d: %d / %d"%(400,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag500 if i>0])
print "%d: %d / %d"%(500,nJumped,nIon)
nJumped = sum([1 for i in jumpedFlag1000 if i>0])
print "%d: %d / %d"%(1000,nJumped,nIon)


pr.prshow()
