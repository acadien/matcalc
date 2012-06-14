#!/usr/bin/python

import sys
import numpy as np
import pylab as pl
#mine
from lmp_cfgIO import readcfg
from paircor import *
from paircorIO import *

def usage():
    print "%s <LAMMPS dump cfg file> <optional: output file>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

fname="paircor_"+sys.argv[1]
if len(sys.argv)==3:
    fname=sys.argv[2]

#Gather the data
bounds,atoms,dummy,types=readcfg(open(sys.argv[1],"r").readlines())
lengths=np.array(map(lambda x: sum(map(fabs,x)),bounds))
atoms=np.add(atoms,np.array(zip(*bounds)[0])*-1) #shift all the atoms over to start at 0,0,0
[rbins,rdist]=paircor_periodic(atoms,lengths,cutoff=10.0,nbins=256)

#Write it to a file
writePairCor(fname,"From file: %s"%sys.argv[1],rbins,rdist)

