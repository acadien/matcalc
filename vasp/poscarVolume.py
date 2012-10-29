#!/usr/bin/python
import sys
#mine
import poscarIO

def unFractional(posfile):
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines(),frac_coord=True)
    for x,y,z in zip(ax,ay,az):
        if x>1.0 or y>1.0 or z>1.0 or x<0.0 or y<0.0 or z<0.0:
            print "Error: it appears these are already unfractional coordinates. %s"%posfile
            exit(0)
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines(),frac_coord=False)
    poscarIO.writeposcar(posfile,[v1,v2,v3], zip(ax,ay,az),atypes,head,frac=False)

def fractional(posfile):
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines(),frac_coord=True)
    fraced=False
    for x,y,z in zip(ax,ay,az):
        if x>1.0 or y>1.0 or z>1.0 or x<0.0 or y<0.0 or z<0.0:
            fraced=True
            break
    if not fraced:
        print "Error: it appears these are already fractional coordiantes. %s"%posfile
        exit(0)
    poscarIO.writeposcar(posfile,[v1,v2,v3], zip(ax,ay,az),atypes,head,frac=True)

def ratio(posfile,ratio):
    print posfile, ratio
    [v1,v2,v3,atypes,ax,ay,az,head,poscar] = poscarIO.readposcar(open(posfile,"r").readlines())
    poscarIO.writeposcar(posfile,[v1,v2,v3],zip(ax,ay,az),atypes,head,True,float(ratio))
    
if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <POSCAR> <frac/unfrac>"%sys.argv[0].split("/")[-1]
        print "%s <POSCAR> <ratio> <ratio value>"%sys.argv[0].split("/")[-1]

    if len(sys.argv)<=2:
        usage()
        exit(0)

    if sys.argv[2]=="frac":
        fractional(sys.argv[1])
    elif sys.argv[2]=="unfrac":
        unFractional(sys.argv[1])
    elif sys.argv[2]=="ratio":
        ratio(sys.argv[1],sys.argv[3])
    else:
        print "Error, unrecognized volume arguement: %s"%sys.argv[2]
        usage()
        exit(0)
