#!/usr/bin/python

import sys
from numpy import dot,cross

def rescaleposcar(poscar,ratio):
    data=open(poscar,"r").readlines()
    ratio=float(ratio)**(1.0/3.0)

    v1=[float(i) for i in data[2].split()]
    v2=[float(i) for i in data[3].split()]
    v3=[float(i) for i in data[4].split()]
    vol=dot(v1,cross(v2,v3))

    print "\n\nUnscaled lattice vectors:"
    print "Volume %g\n"%vol
    for i in data[2:5]:
        print i[0:-2]

    v1p=[i*ratio for i in v1]
    v2p=[i*ratio for i in v2]
    v3p=[i*ratio for i in v3]
    volp=dot(v1p,cross(v2p,v3p))

    data[2]="% 5.18g % 5.18g % 5.18g\n"%(v1p[0],v1p[1],v1p[2])
    data[3]="% 5.18g % 5.18g % 5.18g\n"%(v2p[0],v2p[1],v2p[2])
    data[4]="% 5.18g % 5.18g % 5.18g\n"%(v3p[0],v3p[1],v3p[2])

    print "\n\nScaled Lattice vectors:"
    print "Volume %g\n"%volp
    for i in data[2:5]:
        print i[0:-2]

    val=raw_input("\nHit enter to write back to POSCAR")

    open(poscar,"w").writelines(data)

#The 2nd line of the poscar multiplies each of the lattice constants, 
#altering this value changes the volume
def setScaleMul(poscar,volRatio):
    pcar=open(poscar,"r").readlines()
    ratio=float(volRatio)**(1.0/3.0)
    pcar[1]=str(ratio)+"\n"
    open(poscar,"w").writelines(pcar)

if __name__=="__main__":
    if len(sys.argv)<3:
        print "Usage:"
        print "%s <poscars (space separated)> <volume ratio>"%sys.argv[0]
        exit(0)

    for poscar in sys.argv[1:-1]:
        setScaleMul(poscar,sys.argv[-1])


