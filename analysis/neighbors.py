#!/usr/bin/python
import sys,subprocess
from numpy import *
import pylab as pl
from scipy import weave
from scipy.weave import converters
#mine
from poscar2qvoronoi import poscar2qvoronoi
from voronoiIO import readQVFi,readQVo
from struct_tools import *
from geometry import applyPeriodicBounds
from datatools import flatten

def roundify(atom):
    return map(lambda x:trunc(x*10.01),atom)

def findatom(aval,alist):
    for i,v in enumerate(alist):
        if dist(v,aval)<0.1:
            return i
    raise ValueError("%s not found in list"%str(aval))

#Returns the voronoi neighbor list
def voronoiNeighbors(atoms,basis,atypes=None,style="half"):

    if atypes==None:
        atypes=[0]*len(atoms)
    qvfile,basis,qvatoms,nRealAtoms=poscar2qvoronoi(atoms,basis,atypes)

    #Get hyperplanes & neighbor list from Fi setting of qvoronoi
    p1 = subprocess.Popen(["cat",qvfile],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['qvoronoi','Fi'],stdin=p1.stdout,stdout=subprocess.PIPE).stdout
    qvfidata=p2.readlines()
    polyplanes,allNeighbors=readQVFi(qvfidata)

    #Apply boundaries to ghost atoms
    bounds=array([basis[0][0],basis[1][1],basis[2][2]])
    qvatoms=array(qvatoms)

    #Shift ghosts back into simulation-box and re-form neighborlist
    boundedAtoms=array(map(lambda x:applyPeriodicBounds(x,bounds),qvatoms))
    realAtoms=qvatoms[:nRealAtoms]
    ghost2real=range(nRealAtoms)+[findatom(array(atom),realAtoms) for atom in boundedAtoms[nRealAtoms:]]

    #Build the real neighbor list
    neighbors=[list() for i in range(nRealAtoms)]
    for i,neighbs in enumerate(allNeighbors[:nRealAtoms]):
        ir=ghost2real[i]
        for jNeighb in neighbs:
            jr=ghost2real[jNeighb]
            if jr not in neighbors[ir]:
                neighbors[ir].append(jr)

    if style in ["h","H","half","Half"]:
        return full2half(neighbors)
    return neighbors

#Fast neighbor list generation for orthogonal unit cell vectors
def neighbors(atoms,bounds,r,style="full"):
    #Assumes atom location is >=(min bound) but strictly <(max bound)
    atoms = array(atoms)
    nAtoms = atoms.shape[0]

    stepsz=[0.0,0.0,0.0]
    bounds=[map(float,i) for i in bounds]
    lengths=array([i[1]-i[0] for i in bounds])

    for i,l in enumerate(lengths):
        div=1
        while l/div > r:
            div+=1
            if div>5000:
                print "Error in calculating division size in neighbors() function"
                exit(0)
        stepsz[i]=l/(div-1) if div>1 else l/div

    [Ncx,Ncy,Ncz]=[int(lengths[i]/stepsz[i]) for i in range(3)]
    Ncells=Ncx*Ncy*Ncz
    #Returns the (1-D) index of a cell given its (3-D) coordinates
    def coord2ind((a,b,c)): 
        return int(a*(Ncy*Ncz)+b*Ncz+c)

    #Figure out which atom is in which cell
    def coord2cell(x,y,z):
        return coord2ind(map(int,[ \
                    floor(mod(x,lengths[0])/stepsz[0]),\
                    floor(mod(y,lengths[1])/stepsz[1]),\
                    floor(mod(z,lengths[2])/stepsz[2])]))

    cells=[list() for i in range(Ncells)]
    for i in range(nAtoms):
        cell = coord2cell(atoms[i][0],atoms[i][1],atoms[i][2])
        cells[cell].append(i)
        
    #Returns the set of neighbor cells for celli, includes celli
    def cellNeighbs(celli):
        #Uses Ncell and Ncells, these must be defined correctly as 
        #the total number of cells and cell dimensions respectively.

        coordx=int(floor(celli/(Ncy*Ncz))) #Cell coordinates
        coordy=int(floor(celli%(Ncy*Ncz)/Ncz))
        coordz=int(celli%Ncz) 

        cxs=[(coordx-1)%Ncx, coordx, (coordx+1)%Ncx ] #Neighbor coordinates
        cys=[(coordy-1)%Ncy, coordy, (coordy+1)%Ncy ]
        czs=[(coordz-1)%Ncz, coordz, (coordz+1)%Ncz ]

        cneighbs= [[[(cx,cy,cz) for cz in czs] for cy in cys] for cx in cxs]
        cneighbs= flatten(flatten(cneighbs))
        cneighbs= set([coord2ind(cn) for i,cn in enumerate(cneighbs)])
        return cneighbs

    #Loop over cells, generate neighbor list
    N=len(atoms)
    neighbs=[list() for i in range(N)]
    doneCells=set([ind for ind,cell in enumerate(cells) if len(cell)==0]) #skip all empty cells
    
    for indA in set(range(Ncells)) - doneCells:
        cellA=cells[indA]
        for indB in cellNeighbs(indA) - doneCells:
            cellB=cells[indB]
            for a in cellA:
                
                atomA=atoms[a]
                for b in cellB:
                    atomB=atoms[b]
                    if a==b: continue
                    if dist_periodic(atomA,atomB,lengths) < r:
                        neighbs[a].append(b)
                        if indA!=indB:
                            neighbs[b].append(a)
        doneCells |= set([indA])

    if style in ["h","H","half","Half"]:
        neighbs=full2half(neighbs)
    #print "Done generating Neighbor List."
    return neighbs

#Generate a neighbor list of the N shortest bonds less than length r
def nNearestNeighbors(n,atoms,bounds,r,style="full"):
    neighbs=neighbors(atoms,bounds,r,style)
    bounds=[map(float,i) for i in bounds]
    lengths=array([i[1]-i[0] for i in bounds])

    for i,ineighbs in enumerate(neighbs):
        ds=[[j,dist_periodic(atoms[i],atoms[j],lengths)] for j in ineighbs]
        
        if len(ds)>n:
            ds = sorted(ds,key=lambda x:x[1])[:n]
        
        if len(ds)==0:
            neighbs[i]=[]
        else:
            neighbs[i]=zip(*ds)[0]

    return neighbs

def half2full(hneighbors):
    fneighbors=map(list,hneighbors)
    for i,Neighbs in enumerate(hneighbors):
        for j in Neighbs:
            fneighbors[j].append(i)
    return fneighbors

def full2half(fneighbors):
    hneighbors=map(list,fneighbors)
    for i,a in enumerate(hneighbors):
        for b in a:
            hneighbors[b].remove(i)
    return hneighbors

#Takes a neighbor list (1st shell) and generates the 2nd shell of neighbors for each atom
def secondShell(neighbs):
    shell2=list()

    #loop over each central atoms
    for i,ineighbs in enumerate(neighbs):
        #use set to get uniques
        ashell = list(set(flatten([neighbs[j] for j in ineighbs]+ineighbs)))
        shell2.append(ashell)
    
    return shell2
"""
if __name__=="__main__":
    def usage():
        print "Usage: %s <POSCAR>"%(sys.argv[0])

    if len(sys.argv)<2:
        usage()
        exit(0)

    poscar=open(sys.argv[1],"r").readlines()
    neighbors=voronoiNeighbors(poscar=poscar)

    print neighbors
"""

#This weave code could be optimized to work like the code above.
"""
#neighbBasisCode = ""
int curr=0,ja;
double dist;
for( int i=0; i<na; i++ ){
  for( int j=i+1; j<nall; j++ ){

    dist=0.0;
    for( int x=0; x<3; x++)
      dist += (allatoms[i*3+x]-allatoms[j*3+x]) * (allatoms[i*3+x]-allatoms[j*3+x]);
    dist = sqrt(dist);

    if( dist < rcut ){
      ja=ghostReduce[j];
      neighbFull[i*nbins+(int)neighbCount[i]] = ja;
      neighbFull[ja*nbins+(int)neighbCount[ja]] = i;
      neighbCount[i]++;
      neighbCount[ja]++;
      curr++;
    }
  }
}
return_val = curr;
#""

#Returns the neighbor list for any basis, slower
#Set ghostAtoms=True if the atom list contains ghosts
def neighbors(atoms,basis,rcut,allatoms=None,style='full'):
    print "re-write this to work with weave to speed up!!!!"
    print "if that isn't fast enough, grab the neighbor list from lammps"

    if allatoms==None:
        allatoms,ghostReduce = makeGhosts(atoms,basis)

    nbins=50
    na=len(atoms)
    nall=len(allatoms)
    neighbFull=zeros([nbins*na])-1
    neighbCount=zeros([na])

    nneighbs = weave.inline(neighbBasisCode,['na','nall','nbins','rcut','allatoms','ghostReduce','neighbFull','neighbCount'])
    print na,nall
    neighbFull.shape = [na,nbins]
    neighbFull = [[int(n) for n in set(neighbs) if n>=0] for neighbs in neighbFull]

    if style in ['h','H','half','Half','HALF']:
        neighbs = full2half(neighbFull)
    else:
        neighbs = neighbFull

    #neighbs = [[ghostReduce[j] for j,atomj in enumerate(allatoms) \
    #                if i!=j and dist(atomi,atomj) < rcut ] \
    #               for i,atomi in enumerate(atoms[:500]) ]
    
    return neighbs
"""
