#!/usr/bin/python
import sys
#mine
import poscarIO

def unFractional(posfile):
    [basis,atypes,atoms,head,poscar] = poscarIO.read(open(posfile,"r").readlines(),frac_coord=True)
    for x,y,z in atoms:
        if x>1.0 or y>1.0 or z>1.0 or x<0.0 or y<0.0 or z<0.0:
            print "Error: it appears these are already unfractional coordinates. %s"%posfile
            exit(0)
    [basis,atypes,atoms,head,poscar] = poscarIO.read(open(posfile,"r").readlines(),frac_coord=False)
    poscarIO.write(posfile,basis, atoms,atypes,head,frac=False)

def fractional(posfile):
    [basis,atypes,atoms,head,poscar] = poscarIO.read(open(posfile,"r").readlines(),frac_coord=True)
    fraced=False
    for x,y,z in atoms:
        if x>1.0 or y>1.0 or z>1.0 or x<0.0 or y<0.0 or z<0.0:
            fraced=True
            break
    if not fraced:
        print "Error: it appears these are already fractional coordiantes. %s"%posfile
        exit(0)
    poscarIO.write(posfile,basis,atoms,atypes,head,frac=True)

def ratio(posfile,ratio):
    posdata=open(posfile,"r").readlines()
    posdata[1]="%4.4f\n"%float(ratio)
    open(posfile,"w").writelines(posdata)

if __name__=="__main__":
    def usage():
        print "Usage:"
        print "%s <POSCAR> <frac/unfrac>"%sys.argv[0].split("/")[-1]
        print "%s <POSCAR> <ratio value>"%sys.argv[0].split("/")[-1]

    if len(sys.argv)!=2:
        usage()
        exit(0)

    if sys.argv[2]=="frac":
        fractional(sys.argv[1])
    elif sys.argv[2]=="unfrac":
        unFractional(sys.argv[1])
    else:
        try:
            r=float(sys.argv[2])
            ratio(sys.argv[1],r)
        except:
            print "Error, unrecognized volume arguement: %s"%sys.argv[2]
            usage()
            exit(0)
