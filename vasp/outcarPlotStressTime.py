#!/usr/bin/python
import sys,os
import pylab as pl
#mine
from datatools import windowavg

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <directory that contains OUTCAR file>"

outcar = open(sys.argv[1]+"/OUTCAR","r")

strsx=list()
strsy=list()
strsz=list()

count=0
for line in outcar:
    if "in kB" in line:
        #print line.split("in kB")[1].strip()
        [sx,sy,sz]=[float(i) for i in line.split("in kB")[1].split()[:3]]
        strsx.append(sx)
        strsy.append(sy)
        strsz.append(sz)

#The average Stress
strsavg=[(strsx[i]+strsy[i]+strsz[i])/3 for i in range(len(strsx))]

#A windowed average
winavg=windowavg(strsavg,n=3)

#Plotting
pl.figure()
pl.plot(strsx,ls='--')
pl.plot(strsy,ls='--')
pl.plot(strsz,ls='--')
pl.plot(strsavg)
pl.plot(winavg,c="black",lw=2)
pl.legend(["StressX","StressY","StressZ","Avg","Windowed Avg"],loc=2)
pl.xlabel("Timestep")
pl.ylabel("Stress (kB)")
pl.title(os.getcwd().split("/")[-1])
pl.show()
