#!/usr/bin/python

from numpy import matrix,linalg,asarray,swapaxes
from scipy import dot,cross
import operator
#mine
import poscarIO

def read(elfcar):
    v1,v2,v3,types,xs,ys,zs,header,elfcar=poscarIO.read(elfcar)
    poscardata=(v1,v2,v3,types,xs,ys,zs,header)

    elfcar.pop(0)
    gridsz=[int(i) for i in elfcar.pop(0).split()]

    Tot_pnts = reduce(operator.mul,gridsz)

    elfdata=map(float,"".join(elfcar).split()[:Tot_pnts])
    elfdata=asarray(elfdata)
    elfdata.shape=gridsz
    #elfdata=swapaxes(elfdata,0,2)

    return poscardata,gridsz,elfdata
