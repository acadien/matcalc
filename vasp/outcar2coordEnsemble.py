#!/usr/bin/python

import sys
from math import *
import subprocess
import os
#mine
import outcarIO
import orderParam

from scipy import array,zeros
import pylab as pl

def usage():
    print "Usage: %s <Outcar> "%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)
    
outcarFile = sys.argv[1]
nAtoms = outcarIO.nIons(outcarFile)
basis = array(map(array,outcarIO.basis(outcarFile)))

#Find the starting locations of atomic data in outcarfile
grepResults = subprocess.check_output("grep -b POSITION %s"%outcarFile,shell=True).split("\n")
bytenums=[int(i.split(":")[0]) for i in grepResults if len(i)>2]

outcar= open(outcarFile,"r")
cnOut =open(outcarFile+".cn","w")
cnOut.write("#AverageCN PerAtomCN\n")
for i,b in enumerate(bytenums):
    outcar.seek(b)
    outcar.readline()
    outcar.readline()
    
    atoms = [map(float,outcar.readline().split()[:3]) for a in range(nAtoms)]
    cn,rcut = orderParam.coordinationNumber(atoms,basis)
    cnAvg=sum(cn)/len(cn)
    cns=" ".join(map(str,cn))
    cnOut.write(str(cnAvg)+"\t"+cns+"\n")

#mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4

#opt="vbitrate=2160000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=48 -nosound -o /dev/null mf://\*.png
#mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=2:$opt -mf type=png:fps=48 -nosound -o output.avi mf://\*.png
