#!/usr/bin/python
import sys
#mine
import prshow as pr

import pylab as pl

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0]+" <OUTCAR file>"

outcar = open(sys.argv[1],"r")

tmpt=list()
PE=list()
KE=list()

count=0
for line in outcar:
    if "(temperature" in line:
        tmpt.append(float(line.split("(temperature")[1].split()[0]))
        KE.append(float(line.split("=")[1].split()[0]))
    if "ion-electron   TOTEN" in line:
        PE.append(float(line.split("=")[1].split()[0]))

pl.figure()
pl.plot(tmpt,PE)
pl.plot(tmpt,KE)
pl.plot(tmpt,[PE[i]+KE[i] for i in range(len(PE))])
pl.legend(["PE","KE","TE"])
pl.xlabel("Temperature (K)")
pl.ylabel("Energy (eV)")
pl.title(sys.argv[1])
pr.prshow("outcarET.png")
