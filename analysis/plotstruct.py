import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
#mine
from paircor import paircor, paircor_ang
from duplicate import duplicate26
from struct_tools import *

def plot_atoms(atoms,types,basis=[],dup=0):

	Ntypes=len(set(types))
	N=len(atoms)
	cs=["yellow","green","blue","red","purple","orange","black"]
	
        #Periodic Atoms (duplicated)
	if dup==1:
		atoms,types=duplicate26(atoms,types,zip(v1,v2,v3))

        #Type sorted atoms lists
        tatoms=[list() for i in range(Ntypes)]
	for i,j in enumerate(types):
		tatoms[j-1].append(atoms[i])

	fig=pl.figure()
        aa = fig.add_subplot(111,projection='3d')
        for i in range(Ntypes):
		ax,ay,az=zip(*tatoms[i])
		aa.scatter(ax,ay,az,c=cs.pop(0),marker="o")


	if len(basis)>0:
		v1,v2,v3=zip(*basis)
		aa.plot([0,v1[0]],[0,v1[1]],[0,v1[2]],c="red")
		aa.plot([0,v2[0]],[0,v2[1]],[0,v2[2]],c="blue")
		aa.plot([0,v3[0]],[0,v3[1]],[0,v3[2]],c="green")

def plot_pos_rad_ang(atoms,types,basis,dup=0):

	Ntypes=len(set(types))
	N=len(atoms)
	cs=["yellow","green","blue","red","purple","orange","black"]

        #Periodic Atoms (duplicated, needed for distributions)
        datoms,dtypes=duplicate26(atoms,types,basis)

        #Neighbors needed for angular distribution
        dneighbs=neighbors(datoms,5.0,style="full")

        #Radial Distribution
        [rbins,rdist]=paircor(datoms,inloop=N,cutoff=12.0)

        #Correlation of Angles
        [abins,adist]=paircor_ang(datoms,dneighbs,inloop=N)

        #Type sorted atoms lists
        tatoms=[list() for i in range(Ntypes)]
        if(dup==1):
            for i,j in enumerate(dtypes):
                tatoms[j].append(datoms[i])
        else:
            for i,j in enumerate(types):
                tatoms[j].append(atoms[i])

	v1,v2,v3=zip(*basis)

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
