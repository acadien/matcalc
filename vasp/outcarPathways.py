#!/usr/bin/python
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import sys
from math import *
from scipy import array,zeros
import pylab as pl
from scipy import weave
from scipy.weave import converters
#mine
import outcarIO
import meanSquareDist
from struct_tools import dist
#Calculates mean squared distance

def usage():
    print "Usage: %s <Outcar> <Config0> "%sys.argv[0].split("/")[-1]

def outcarPaths(outcarfile):
    outcar=open(outcarfile,"r")
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

    atoms=array(atoms).ravel()
    delT=array(range(1,Ntime+1))
    compiler_args=['-march=native -O3 -fopenmp']
    headers=r"""#include <omp.h>"""
    libs=['gomp']
    weave.inline(meanSquareDist.undoPBCcode,['atoms','Ntime','Natom','lengths'],\
                     extra_compile_args=compiler_args,\
                     support_code=headers,\
                     libraries=libs)

    
    distances=[dist(atoms[0,i],atoms[-1,i]) for i in range(Natom)]
    print distances 
    exit(0)
    atoms.shape=(Ntime,Natom,3)

    fig=pl.figure()
    ax=fig.add_subplot(111,projection='3d')
    for i in range(Natom):
        ax.plot(atoms[:,i,0],atoms[:,i,1],atoms[:,i,2])

    pl.suptitle(outcarfile)
    pl.show()

if __name__=="__main__":
    if len(sys.argv)<2:
        usage()
        exit(0)
    
    outcarfile=sys.argv[1]
    try:
        num=sys.argv[1].split("_")[1]
    except IndexError:
        num=None

    outcarPaths(outcarfile)


