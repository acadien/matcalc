#!/usr/bin/python
import sys
import pylab as pl

if len(sys.argv)<2:
    print "Usage:"
    print sys.argv[0].split("/")[-1]+" <OUTCAR file>"

outcar = open(sys.argv[1],"r")

PE=list()
KE=list()

count=0
for line in outcar:
    if "(temperature" in line:
        KE.append(float(line.split("=")[1].split()[0]))
    if "ion-electron   TOTEN" in line:
        PE.append(float(line.split("=")[1].split()[0]))

N=len(PE)
pl.figure()
pl.plot(PE)
pl.plot([PE[i]+KE[i] for i in range(len(PE))])
pl.legend(["PE","TE"])
pl.xlabel("Timestep")
pl.ylabel("Energy (eV)")
pl.title(sys.argv[1])
pl.show()
