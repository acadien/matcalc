#!/usr/bin/python
import pylab as pl
from numpy import *
import sys

def usage():
    print "usage: %s <run number> "%sys.argv[0].split("/")[-1]

#expects 9 files:
#6 vasp files (2 at each pressure: energy & stress)
#3 lammps dump files (1 at each pressure)

if len(sys.argv)<2:
    usage()
    exit(0)

#file information
run=int(sys.argv[1])
natom=216
vaspfileE="feedback%d_energy"%run
vaspfileS="feedback%d_stress"%run
lmpssuffix="GPa_nvt%d_thermo.dat"%run
pressures=["00","20","40"]

#Parse the VASP data
def parseVaspE(line):
    line=line.split()
    fname=line[0].split("/")[0]
    toten=float(line[-2])/natom
    pressure=fname.split("_")[0].strip("GPa")
    num=int(fname.split("_")[-1])
    cool=0
    if fname.split("_")[1].strip("0")=="cool": cool=1
    return [pressure,cool,num,toten]

def parseVaspS(line):
    line=line.split()
    fname=line[0].split("/")[0]
    stress=float(line[-7])/10.
    pressure=fname.split("_")[0].strip("GPa")
    num=int(fname.split("_")[-1])
    cool=0
    if fname.split("_")[1].strip("0")=="cool": cool=1
    return [pressure,cool,num,stress]

vPress,vCool,vNum,vStress=map(list,zip(*sorted(map(parseVaspS,open(vaspfileS,"r").readlines()))))
vPress,vCool,vNum,vToten=map(list,zip(*sorted(map(parseVaspE,open(vaspfileE,"r").readlines()))))

#Convert numbering and heat/cool info to recoup the iteration number from lammps simulation
vIter=[int(c*1E6+2E5+n*5E4) for c,n in zip(vCool,vNum)]

#Split up total energy and iterations by pressure changes
s2= [i+1 for i,[a,b] in enumerate(zip(vPress,vPress[1:]+[vPress[0]])) if a!=b][:-1]
vToten=[vToten[:s2[0]],vToten[s2[0]:s2[1]],vToten[s2[1]:]]
vStress=[vStress[:s2[0]],vStress[s2[0]:s2[1]],vStress[s2[1]:]]
vIter=[vIter[:s2[0]],vIter[s2[0]:s2[1]],vIter[s2[1]:]]

#Parse the LAMMPS data, grab only iterations from vasp
lToten=list()
lTemp=list()
lPress=list()

for p,iters in zip(pressures,vIter):
    iters=map(str,iters)
    temp,tE,pres=list(),list(),list()
    rdy=False
    for line in open(p+lmpssuffix,"r").readlines():
        if "Step TotEng Temp Press Volume" in line:
            rdy=True
        if not rdy:
            continue
        try:
            line=line.split()
        except IndexError:
            continue
        if len(line)>0 and line[0] in iters:
            iters.remove(line[0])
            temp.append(line[2])
            tE.append(float(line[1])/natom)
            pres.append(float(line[3])/10000.)
    lToten.append(tE)
    lTemp.append(temp)
    lPress.append(pres)

i=0
lTemp=map(lambda x:map(lambda y:round(float(y)/10,0)*10,x),lTemp)
#print lTemp
pl.figure()
for p,vi,ve,vs,le,lt,ls in zip(pressures,vIter,vToten,vStress,lToten,lTemp,lPress):
    i+=1
    pl.subplot(320+i)
    if i<=1:
        pl.title("                                                        Feedback %d"%run)
    pl.plot(ve,label=p+"GPa VASP")
    pl.plot(le,label="LAMMPS")
    pl.ylabel("eV/atom")
    pl.legend()
    #print array(ve)-array(le)
    #print sum(array(ve[1:])-array(le[1:]))/len(le[1:])
    #print sum(array(vs[1:])-array(ls[1:]))/len(ls[1:])
    i+=1
    pl.subplot(320+i)
    pl.plot(vs,label=p+"GPa")
    pl.plot(ls,label="LAMMPS")
    pl.ylabel("Stress GPa")
    pl.legend()
    
pl.show()

