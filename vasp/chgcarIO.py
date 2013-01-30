#!/usr/bin/python

import scipy
from scipy import *
import operator
#mine
import poscarIO


#poscardata,gridsz,chgdata = chgcarIO.read(open(chgcarFileName,"r").readlines())
#sum(sum(sum(chgdata))) = Total # of electrons (or valence electrons)
def read(chgcar,frac_coord=True):
    basis,types,atoms,header,chgcar=poscarIO.read(chgcar,frac_coord)
    poscardata=(basis,types,atoms,header)

    chgcar.pop(0)
    gridsz=[int(i) for i in chgcar.pop(0).split()]

    Tot_pnts = reduce(operator.mul,gridsz)

    #Divices by Tot_pnts to get absolute Charge (e.g. total # electrons)
    chgdata=array([float(i)/Tot_pnts for i in "".join(chgcar).split()[:Tot_pnts]])
    chgdata=asarray(chgdata)
    chgdata.shape=gridsz
    chgdata=swapaxes(chgdata,0,2)

    return poscardata,gridsz,chgdata

#Writes a chargecar to vtk format
def writeVTK(fname,poscardata,gridsz,chgdata,clean=False):
    basis,types,atoms,header=poscardata
    
    Tot_pnts=reduce(operator.mul,gridsz)

    header=list()
    header.append("# vtk DataFile Version 2.0")
    header.append("Charge Density 3D grid")
    header.append("ASCII")
    header.append("DATASET STRUCTURED_POINTS")
    header.append("DIMENSIONS "+" ".join(map(str,gridsz)))
    header.append("ORIGIN 0.0 0.0 0.0")
    header.append("SPACING "+" ".join(map(str,[i/j for i,j in zip([basis[0][0],basis[1][1],basis[2][2]],gridsz)])))
    header.append("POINT_DATA %d"%Tot_pnts)
    header.append("SCALARS charge_density float")
    header.append("LOOKUP_TABLE default")
    
    vtkfile=open(fname,"w")
    vtkfile.writelines([x+"\n" for x in header])

    if clean==True:
        chgdata=chgdata.ravel()
        def clearoutlier(i,a): return 1.0 if i>a else i/a
        avgval=scipy.sum(chgdata)/Tot_pnts
        vclearout=vectorize(clearoutlier)
        chgdata=vclearout(chgdata,avgval)
        #chgdata=array([i/avgval for i in chgdata])
        chgdata.shape=gridsz
    #chgdata=chgdata.ravel()
    #savetxt(fname,chgdata,delimiter=' ')
    for i in range(gridsz[0]):
         for j in range(gridsz[1]):
            vtkfile.write(" ".join(map(lambda x:str(x),chgdata[i][j]))+"\n")
