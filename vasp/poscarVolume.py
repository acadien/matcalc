#!/usr/bin/python
import sys
from numpy import matrix,linalg
#mine
poscar

def unFractional(posfile):
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines())
    poscarIO.writeposcar(posfile,[v1,v2,v3], zip(ax,ay,az),atypes,head,frac=False)

def fractional(posfile):
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines(),frac_coord=True)
poscarIO.writeposcar(sys.argv[1],[v1,v2,v3], zip(ax,ay,az),atypes,head,frac=False)

def ratio(posfile):
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines())
    poscarIO.writeposcar(posfile,[v1,v2,v3], zip(ax,ay,az),atypes,head,frac=True,ratio)
    
if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <POSCAR> <frac/unfrac/ratio> <ratio value>"%sys.argv[0].split("/")[-1]

    if len(sys.argv)<=2:
        usage()
        exit(0)

    if sys.argv[2]=="frac":
        fractional(sys.argv[1])
    elif sys.argv[2]=="unfrac":
        unFractional(sys.argv[1])
    elif sys.argv[2]=="ratio"
        ratio(sys.argv[1],sys.argv[2])
    else:
        print "Error, unrecognized fraction request"
        usage()
        exit(0)
