#!/usr/bin/python
import sys
from math import *
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
import os
#mine
#from plotstruct import plot_pos_rad_ang
from rdf import rdf,adf
from duplicate import duplicate26
from struct_tools import *
import poscarIO
from datatools import flatten

def findsymmetry(vs,ts,xyzpnts):
    tfname="tmpsym.in"
    tf=open(tfname,"w")
    wrdata="Comment line\n"
    wrdata+="0\taccuracy\n"
    wrdata+="1\tlattice form: vectors\n"
    wrdata+="%g %g %g\n"%(vs[0][0],vs[0][1],vs[0][2])
    wrdata+="%g %g %g\n"%(vs[1][0],vs[1][1],vs[1][2])
    wrdata+="%g %g %g\n"%(vs[2][0],vs[2][1],vs[2][2])
    wrdata+="1\n"
    wrdata+="1 0 0\n0 1 0\n0 0 1\n"
    wrdata+="%d\n"%(sum(ts))
    types=" ".join([str(k) for k in flatten([[i+1 for j in range(ts[i])] for i in range(len(ts))])])
    wrdata+="%s\n"%types
    for i,dat in enumerate(xyzpnts):
        wrdata+="%g %g %g\n"%(dat[0],dat[1],dat[2])
    tf.write(wrdata)
    tf.close()

    #Issue the findsym command and parse the data
    os.system("findsym < %s > temp"%tfname)
    os.system("rm %s"%tfname)
    utpt=open("temp","r").readlines()
    os.system("rm temp")
    for line in utpt:
        if "space_group_name" in line:
            print "Symmetry Space Group= "+" ".join(line.split()[1:])
            break


#########################################
#                 MAIN                  #
#########################################

#Reading in the options and preparing for file reading
if len(sys.argv) != 4:
    print "Usage: %s <start_gen=1> <plot_dup=0>"%sys.argv[0]
    print "Must be run from the USPEX results directory.\n"
elif len(sys.argv) != 3:
    print "Usage: %s <start_gen=1> <plot_dup=0>"%sys.argv[0]
    print "Must be run from the USPEX results directory."
    exit(0)

start=1
if len(sys.argv) >= 2:
    start=int(sys.argv[1])

pdup=0
if len(sys.argv) == 3:
    pdup=int(sys.argv[2])

colors=["red","blue","green","yellow","orange","purple","black"]
poscar = open("./BESTgatheredPOSCARS","r").readlines()
fitness = open("./BESTfitness","r").readlines()
kpoints = open("./BESTkpoints.dat","r").readlines()
volumes = open("./BESTvolumes.dat","r").readlines()
enthalpy = open("./BESTenthalpies.dat","r").readlines()
gencnt=0

colors=["red","blue","green","yellow","orange","purple","black"]
while True:
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.read(poscar)

    if v1==v2==v3==-1:
        break

    gencnt+=1
    
    eanum=int(head.split()[0].split("EA")[1])

    j=0
    types=list()
    cs=list()
    for i in atypes:
        types+=[j]*i
        cs+=[colors[j]]
        j+=1
    N=len(types)
    Ntypes=len(set(types))

    print "="*80
    print "Generation %d, EA#%d"%(gencnt,eanum)
    a=fitness.pop(0).strip()
    print "Fitness= %s eV   PerAtom Energy= %s eV\nkpoints= %sVolume= %sEnthalpy= %s"%(a,round(float(a)/N,3),kpoints.pop(0),volumes.pop(0),enthalpy.pop(0).strip())

    if gencnt >= start:
        #Actual Atoms
        atoms=zip(ax,ay,az)

        #plot_pos_rad_ang(atoms,types,pdup=pdup)

        #Periodic Atoms (duplicated)
        datoms,dtypes,dbasis=duplicate26(atoms,types,zip(v1,v2,v3))

        #Neighbors needed for angular distribution
        dneighbs=neighbors(datoms,5.0,style="full")

        #Radial Distribution
        [rbins,rdist]=rdf(datoms,inloop=N,cutoff=12.0)

        #Correlation of Angles
        [abins,adist]=adf(datoms,dneighbs,inloop=N)

        #Type sorted atoms lists
        tatoms=[list() for i in range(Ntypes)]
        if(pdup==1):
            for i,j in enumerate(dtypes):
                tatoms[j].append(datoms[i])
        else:
            for i,j in enumerate(types):
                tatoms[j].append(atoms[i])

        #findsymmetry([v1,v2,v3],types,zip(xs,ys,zs))

        fig=pl.figure(figsize=(12,6))
        aa = fig.add_subplot(311,projection='3d')
        aa.set_position([0,0,0.5,1.0])
        for i in range(Ntypes):
            ax,ay,az=zip(*tatoms[i])
            aa.scatter(ax,ay,az,c=cs.pop(0),marker="o")

        aa.plot([0,v1[0]],[0,v1[1]],[0,v1[2]])
        aa.plot([0,v2[0]],[0,v2[1]],[0,v2[2]])
        aa.plot([0,v3[0]],[0,v3[1]],[0,v3[2]])
        
        bb = fig.add_subplot(312)
        bb.set_position([0.56,0.55,0.4,0.4])
        bb.plot(rbins,rdist)
        pl.xlabel("Radius (A)")
        pl.ylabel("Count")

        cc = fig.add_subplot(313)
        cc.set_position([0.56,0.08,0.4,0.36])
        cc.plot(abins,adist)
        pl.xlabel("Angle (deg)")
        pl.ylabel("Count")

        pl.show()

