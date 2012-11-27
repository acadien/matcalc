#!/usr/bin/python
import sys
#mine
import outcarIO

def outcar2poscar(outcar,outposcar,wantconfig):
    pcar=outcarIO.outcar2poscar(outcar,wantconfig)
    open(outposcar,"w").writelines(pcar)

if __name__ == "__main__":
    if len(sys.argv)<3:
        print "Usage:"
        print sys.argv[0]+" <OUTCAR> <output poscar name> <desired config #, 0-indexed (default is last)>"
        exit(0)

    outcar=sys.argv[1]
    outposcar=sys.argv[2]

    if len(sys.argv)==4:
        wantconfig=int(sys.argv[3])
    else:
        wantconfig=-1 #last one by default
        
    outcar2poscar(outcar,outposcar,wantconfig)
