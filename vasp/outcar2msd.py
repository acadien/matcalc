#!/usr/bin/python
import sys
from math import *
from scipy import array,zeros
import pylab as pl
#mine
from meanSquareDist import meanSquareDistRef
import outcarIO
#msdIO

#Calculates mean squared distance

def usage():
    print "Usage: %s <Outcar> <Config0> "%sys.argv[0].split("/")[-1]

def outcarMeanSquareDisplace(outcarFile,refStructure=None):
    outcar=open(outcarFile,"r")
    atoms=list() #atoms[time][atom]

    #Grab ion types
    while True:
        line=outcar.readline()
        if "ions per type" in line:
            break
    atypes=map(int,line.split("=")[1].split())
    Natom=sum(atypes)

    #Grab basis vectors
    while True:
        line=outcar.readline()
        if "direct lattice vectors" in line:
            break
    basis=array([map(float,outcar.readline().split()[:3]) for i in range(3)])
    lengths=array([basis[0][0],basis[1][1],basis[2][2]])

    #Grab atom positions
    count=0
    posit=False
    for line in outcar:
        if posit:
            if "--" in line:
                if len(a)==0:
                    continue
                else:
                    #Analysis
                    atoms.append(array(a))
                    posit=False
            else:
                a.append(map(float,line.split()[:3]))
        elif "POSITION" in line:
            a=list()
            posit=True
            count+=1
    atoms=array(atoms)
    Ntime=len(atoms)

    if refStructure==None:
        delT,msd=meanSquareDistRef(atoms,0,Natom,Ntime,lengths)
    elif refStructure > Ntime:
        print "%d requested but %d structures max"%(refStructure,Ntime)
    else:
        delT,msd=meanSquareDistRef(atoms,refStructure,Natom,Ntime,lengths)

    return delT,msd

if __name__=="__main__":
    if len(sys.argv)<2:
        usage()
        exit(0)
    
    outcarFile=sys.argv[1]
    try:
        num=sys.argv[1].split("_")[1]
    except IndexError:
        num=None

    refStructure=None
    if len(sys.argv)==3:
        refStructure=int(sys.argv[2])

    delT,msd=outcarMeanSquareDisplace(outcarFile,refStructure)

    if refStructure==None:
        msdfile=outcarFile+".msd"
    else:
        msdfile=outcarFile+".msd%d"%refStructure

    header=["Mean Squared Displacement from %s.\n"%(sys.argv[0]),"delTimeStep msd($\AA^2$)\n"]
    msddata=header+[str(x)+" "+str(y)+"\n" for x,y in zip(delT,msd)]
    print "Writing %s."%msdfile
    open(msdfile,"w").writelines(msddata)

#    pl.plot(msd)
#    pl.show()

