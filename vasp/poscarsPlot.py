#!/usr/bin/python
import sys
from math import *
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt
#mine
import poscarIO
from duplicate import duplicate26

def plotsimulation(basis,atoms,types,fig=None):
    colors=["red","blue","green","yellow","orange","purple","black"]
    if fig==None:
        fig=plt.figure()
    aa = fig.gca(projection='3d')
    ax,ay,az=zip(*atoms)
    ind=0
    for cindex,i in enumerate(types):
        aa.scatter(ax[ind:ind+i],ay[ind:ind+i],az[ind:ind+i],c=colors[cindex],marker="o")
        ind+=i
    for i in basis:
        aa.plot([0,i[0]],[0,i[1]],[0,i[2]])
    return aa

if __name__=="__main__":

    #Reading in the options and preparing for file reading
    if len(sys.argv)<2:
        print "Usage:"
        print sys.argv[0]+" <file with a bunch of poscars> <duplicate=0>"
        print "If duplicate==1, copies the structure (translated by lattice vecs) and then plots."
        exit(0)


    poscar = open(sys.argv[1],"r").readlines()

    dupl=0
    if len(sys.argv)==3:
        dupl=int(sys.argv[2])

    while True:
        if len(poscar)<=1:
            break
        [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.read(poscar)
        types=list()
        #cs=list()
        j=1
        for i in atypes:
            types+=[j]*i
            #cs.append(colors[j-1])
            j+=1
        Ntypes=len(types)

        atoms=zip(ax,ay,az)

        basis=[v1,v2,v3]

        if dupl==1:
            atoms,types,basis = duplicate26(atoms,types,zip(v1,v2,v3))
            atoms,types = zip(*sorted(zip(atoms,types),key=lambda x:x[1])) #sort based on type
            atypes = [i*26 for i in atypes]           
            #ax,ay,az=zip(*atoms)          
            #basis=[[j*3 for j in i] for i in basis]
            print basis

        plotsimulation(basis,atoms,atypes)
        pl.show()
