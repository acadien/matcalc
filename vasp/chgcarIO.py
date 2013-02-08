#!/usr/bin/python

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
    chgdata=array(map(lambda x:float(x)/Tot_pnts,"".join(chgcar).split()))

    #To eliminate points that step outside a reasonable range
    avg = sum(chgdata)/len(chgdata)

    #chgdata[where(chgdata > 3*avg)[0]]=avg

    chgdata.shape=gridsz
    #Some chgcars seem to need their axes swapped, possible bug in VASP
#    chgdata=swapaxes(chgdata,0,2)

    return poscardata,chgdata

#Writes a chargecar to vtk format
def writeVTK(fname,poscardata,chgdata):
    basis,types,atoms,header=poscardata
    gridsz=chgdata.shape

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

    mx=chgdata.max()
#    mx=1
    for i in range(gridsz[0]):
         for j in range(gridsz[1]):
            vtkfile.write(" ".join(map(lambda x:str(x/mx),chgdata[i][j]))+"\n")
