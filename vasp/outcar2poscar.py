#!/usr/bin/python
import sys
#mine
import outcarIO
import poscarIO

def outcar2poscar(outcar,outposcar,wantconfig):
    basis,atypes,atoms,head,poscar = outcarIO.outcar2poscar(outcar,wantconfig)
    poscarIO.write(outposcar,basis,atoms,atypes,head)

if __name__ == "__main__":
    if len(sys.argv)<3:
        print "Usage:"
        print sys.argv[0].split("/")[-1]+" <OUTCAR> <output POSCAR> <desired config #, 0-indexed (default is last)>"
        exit(0)

    outcar=sys.argv[1]
    outposcar=sys.argv[2]

    if len(sys.argv)==4:
        wantconfig=int(sys.argv[3])
    else:
        wantconfig=-1 #last one by default
        
    outcar2poscar(outcar,outposcar,wantconfig)
