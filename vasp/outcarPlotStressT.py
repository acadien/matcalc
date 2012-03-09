#!/usr/bin/python
import sys
import pylab as pl
#mine
from datatools import windowavg

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <directory that contains OUTCAR file>"

outcar = open(sys.argv[1]+"/OUTCAR","r")

tmpt=list()
strsx=list()
strsy=list()
strsz=list()

count=0
for line in outcar:
    if "(temperature" in line:
        tmpt.append(float(line.split("(temperature")[1].split()[0]))
    if "in kB" in line:
        #print line.split("in kB")[1].strip()
        [sx,sy,sz]=[float(i) for i in line.split("in kB")[1].split()[:3]]
        strsx.append(sx)
        strsy.append(sy)
        strsz.append(sz)

#The average Stress
strsavg=[(strsx[i]+strsy[i]+strsz[i])/3 for i in range(len(strsx))]

#A windowed average
winavg=windowavg(strsavg,10)

print len(winavg),len(strsavg),len(strsz),len(tmpt)

#Plotting
pl.figure()
pl.plot(tmpt,strsx,ls='--')
pl.plot(tmpt,strsy,ls='--')
pl.plot(tmpt,strsz,ls='--')
pl.plot(tmpt,strsavg,ls='--')
pl.plot(tmpt,winavg,c="black",lw=2)
pl.legend(["StressX","StressY","StressZ","Avg","Running Avg"],loc=4)
pl.xlabel("Temperature (K)")
pl.ylabel("Stress (kB)")
pl.title(sys.argv[1])
pl.show()
