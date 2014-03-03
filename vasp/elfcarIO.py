#!/usr/bin/python

from numpy import matrix,linalg,asarray,swapaxes
import numpy as np
from scipy import dot,cross
import operator
#mine
import poscarIO

def read(elfcar):
    basis,types,atoms,header,elfcar=poscarIO.read(elfcar)

    atoms=asarray(atoms)
    basis=asarray(basis)
    poscardata=(basis,types,atoms,header)

    elfcar.pop(0)
    gridsz=[int(i) for i in elfcar.pop(0).split()]
    gridsz=[gridsz[2],gridsz[1],gridsz[0]]
    Tot_pnts = reduce(operator.mul,gridsz)

    elfdata=asarray(map(float,"".join(elfcar).split()[:Tot_pnts]))
    elfdata.shape=gridsz
    elfdata=elfdata.swapaxes(0,2)
    return poscardata,elfdata


