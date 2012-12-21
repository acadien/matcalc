#!/usr/bin/python

import sys,os,time
from math import degrees
from numpy import matrix,linalg
#mine
import poscarIO
from struct_tools import vecs2lattice

def usage():
    print "%s <element1> <element2> ... <elementN> <in:POSCAR-file> <out:CIF-file>"%(sys.argv[0])

if len(sys.argv)<4:
    usage()
    exit(0)

elems=sys.argv[1:-2]
try:
    poscar=open(sys.argv[-2],"r").readlines()
    xyz=open(sys.argv[-1],"w")
except IOError:
    print "Error: Unable to open files."
    usage()
    exit(0)

[basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)

N=sum(atypes)

#Calculate fractional coordinates of atoms
A=matrix(basis)
[fx,fy,fz]=[list(i) for i in zip(*[linalg.solve(A,p)[:] for p in atoms])]

#Calculate lattice parameters
a,b,c,A,B,C=vecs2lattice(*basis)
A,B,C=degrees(A),degrees(B),degrees(C)

data= "_audit_creation_date              %s\n"%time.ctime()
data+="_pd_phase_name                    \'poscar2cif.py from %s/%s\'\n"%(os.getcwd(),sys.argv[-2])
data+="_cell_length_a                    %5f\n"%a
data+="_cell_length_b                    %5f\n"%b
data+="_cell_length_c                    %5f\n"%c
data+="_cell_angle_alpha                 %5f\n"%A
data+="_cell_angle_beta                  %5f\n"%B
data+="_cell_angle_gamma                 %5f\n"%C
data+="_symmetry_space_group_name_H-M    \'P 1\'\n"
data+="_symmetry_Int_Tables_number       1\n"
data+="loop_\n"
data+="_symmetry_equiv_pos_as_xyz\n"
data+="    \'x, y, z\'\n"
data+="loop_\n"
data+="   _atom_site_type_symbol\n"
data+="   _atom_site_label\n"
data+="   _atom_site_fract_x\n"
data+="   _atom_site_fract_y\n"
data+="   _atom_site_fract_z\n"

for ntype in range(len(atypes)):
    cnt=0
    for i in range(atypes[ntype]):
        cnt+=1
        data+="   %s   %s%3.3d\t%10.5f\t%10.5f\t%10.5f\n"%(elems[ntype],elems[ntype],cnt,fx.pop(0),fy.pop(0),fz.pop(0))

xyz.write(data)
