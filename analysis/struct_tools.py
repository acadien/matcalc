#!/usr/bin/python

#A collection of helper functions and general utilities for analzing a set of atoms

from math import *
from scipy import weave
from scipy.weave import converters
from numpy import *
import pylab as pl
#mine
from datatools import flatten       

#===========  Helper Geometry Functions  ================
#returns the scalar distance between the 2 lists/arrays p1 and p2

distcode = """
double c=0.0;
for(int i=0;i<3;i++)
  c = c + (a[i]-b[i])*(a[i]-b[i]);
return_val=sqrt(c);
"""
def dist(a,b):
    return weave.inline(distcode,['a','b'])

distcodeP = """
double c=0.0,d=0.0;x
for(int i=0;i<3;i++){
  d = a[i]-b[i];
  if(d>lengths[i]/2.0) d -= lengths[i];
  if(d<-lengths[i]/2.0) d += lengths[i];
  c = c + d*d;
}
return_val=sqrt(c);
"""
#Calculates dist between 2 points enforcing periodic boundary conditions
def dist_periodic(a,b,lengths):
    return weave.inline(distcodeP,['a','b','lengths'])

distSPcode = """
double d=0.0,c=0.0;
for(int tx=-1;tx<2;tx++){
    for(int ty=-1;ty<2;ty++){
        for(int tz=-1;tz<2;tz++){
            c=0.0;
            d=a[0]-b[0]+tx;
            c+=d*d;
            d=a[1]-b[1]+ty;
            c+=d*d;
            d=a[2]-b[2]+tz;
            c+=d*d;
            rs[(tx+1)*9+(ty+1)*3+tz+1]=c;
        }
    }
}
"""
def distsp(a,b,l,rs):
    weave.inline(distSPcode,['a','b','l','rs'])
    return rs

angcode = """
double ang,x,dab=0.0,dac=0.0,dbc=0.0;
for(int i=0;i<3;i++)
  dab = dab + (a[i]-b[i])*(a[i]-b[i]);
for(int i=0;i<3;i++)
  dac = dac + (a[i]-c[i])*(a[i]-c[i]);
for(int i=0;i<3;i++)
  dbc = dbc + (b[i]-c[i])*(b[i]-c[i]);
x=(dab + dbc - dac)/(2.0*sqrt(dab)*sqrt(dbc));
if(fabs(fabs(x)-1.0) <= 1e-9)
  return_val=0.0;
else
  return_val=acos(x);
"""
def ang(a,b,c):
    return weave.inline(angcode,['a','b','c'])

def mag(vec):
    return sqrt(dot(vec,vec))

def normalize(vec):
    a=mag(vec)
    return vec/a

def rotmatx(vec1,vec2):#Rotates vector1 into vector2 and returns the rotation matrix
    raxis=normalize(cross(vec1,vec2))
    [rx,ry,rz]=raxis
    rtheta=ang(vec1,array([0,0,0]),vec2)
    cth=cos(rtheta)
    sth=sin(rtheta)

    RM=array([\
        [cth+rx*rx*(1-cth),rx*ry*(1-cth)-rz*sth,rx*rz*(1-cth)+ry*sth],\
        [ry*rx*(1-cth)+rz*sth,cth+ry*ry*(1-cth),ry*rz*(1-cth)-rx*sth],\
        [rz*rx*(1-cth)-ry*sth,rz*ry*(1-cth)+rx*sth,cth+rz*rz*(1-cth)]\
        ])

    return RM

def init():
    #atoms=map(array,[[0.5, 1.5, 2.5], [4.9, 4, 2], [1.1, 1.5, 2.5]])
    atoms=map(array,[[0.5,0.5,0.5], [1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [0.5,0.5,1.0],[0.5,1.0,1.0],[1.0,0.5,1.0],[1.0,1.0,0.5],[1.0,1.0,1.0]])
    bounds=[[0,5],[0,5],[0,5]]
    r=0.99
    return atoms,bounds,r

def volume(a,b,c):
    return dot(a,cross(b,c))

#===========  Neighbor List  ===================
def neighbors(atoms,bounds,r,style="full"):
    #Assumes atom location is >=(min bound) but strictly <(max bound)

    #Makes a neighbor list for the list of atoms
    #Neighbors within a distance r are considered
    #Added functionality for including periodicity

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
    coord2ind=lambda (a,b,c): int(a*(Ncy*Ncz)+b*Ncz+c)

    #Figure out which atom is in which cell
    coord2cell=lambda at: coord2ind(map(int,[floor(at[0]/stepsz[0]),floor(at[1]/stepsz[1]),floor(at[2]/stepsz[2])]))
    
    cells=[list() for i in range(Ncells)]
    for ind,cell in enumerate(map(coord2cell,atoms)):
        cells[cell].append(ind)
        
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
            #print indA,indB
            #print "="*50
            cellB=cells[indB]

            #Loop over atoms in each cell
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
        for i,a in enumerate(neighbs):
            for b in a:
                neighbs[b].remove(i)
    return neighbs


def old_neighbors(atoms,bounds,r,style="full"):
    lengths=array(map(lambda x: x[1]-x[0],bounds))

    N=len(atoms)
    neighbs=[list() for i in range(N)]
    if style in ["h","H","half","Half"]:
        neighbs = [ [j for j in range(i+1,N) if dist_periodic(atoms[i],atoms[j],lengths) < r] for i in range(N)]
    elif style in ["f","F","Full","full"]:
        for i,ai in enumerate(atoms):
            for j,aj in enumerate(atoms[i+1:]):
                j+=i+1
                if dist_periodic(ai,aj,lengths) < r:
                    neighbs[i].append(j)
                    neighbs[j].append(i)
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


#========  Convert Cell Vector basis  ==========
def vecs2lattice(v1,v2,v3):
    #Given a 3 vector basis v1,v2,v3
    #Returns the lattice lengths (a,b,c) and angles (A,B,C)
    v1=array(v1)
    v2=array(v2)
    v3=array(v3)
    origin=array([0,0,0])

    a=(v1[0]**2+v1[1]**2+v1[2]**2)**0.5
    b=(v1[0]**2+v1[1]**2+v1[2]**2)**0.5
    c=(v1[0]**2+v1[1]**2+v1[2]**2)**0.5

    A=ang(v1,origin,v2)
    B=ang(v1,origin,v3)
    C=ang(v2,origin,v3)
    return a,b,c,A,B,C


#===========  Find Structures ==================
def find_squares(atoms,neighbs,atol=1.0,dtol=0.1):#tolerance in degrees
    atol=radians(atol)
    N=len(atoms)
    right=pi/2.0
    squares=list()
    uniqs=list()

    for i,ai in enumerate(atoms):
        for j in neighbs[i]:
            aj=atoms[j]
            dij=dist(ai,aj)
            for k in neighbs[i]:
                if j==k: continue
                ak=atoms[k]
                angle = ang(ai,aj,ak)
                if fabs(dist(aj,ak)-dij) > dtol: continue
                if fabs(angle - right) <= atol:
                    fourth = set(neighbs[k]) - set([i,j,k])
                    for l in fourth:
                        al = atoms[l]
                        if fabs(dist(ai,al)-dij) > dtol: continue
                        if fabs(dist(ak,al)-dij) > dtol: continue
                        square=[i,j,k,l]
                        ssq = sorted(square)
                        if ssq not in uniqs:
                            squares.append(square)
                            uniqs.append(ssq)

    return squares
     
def find_cubes(squares):
    multiverts=[v for v in flatten(squares)]
    vertices=[v for v in set(multiverts)]

    #These vertices are a part of 3 different squares 
    rightverts = [v for v in vertices if multiverts.count(v) >= 3]

    #Check if each value is in some biglist
    def allin(values,biglist):
        for value in values:
            if value not in biglist:
                return False
        return True
    
#    for i,isq in enumerate(squares):
#        if allin(isq,rightverts):
#            print isq


#================= Making Polyhedron from points (Voronoi) ============

#Given a list of points all in the same hyperplane, defined by the coefs from 'hyperplane' return the list of points sorted and ready to be used as a polygon.
def points2polygon(points,hyperplane):
    com=array([sum(i)/len(i) for i in zip(*points)])
    points2=[i-com for i in points]
    s=min([i for i in zip(range(len(points2)),points2) if i[1][2]>0],key=lambda x:x[1][2])
    rot=rotmatx(points2[0],array([1,0,0]))
    points2=[dot(rot,p) for p in points2]
    com=zeros(3)
    start=points2[0]
    angles=[ang(start,com,i) for i in points2]
    for i,p in enumerate(points2):
        if p[2]<0 and points2[0][2]>0:
            angles[i]=2*pi-angles[i]
        if p[2]>0 and points2[0][2]<0:
            angles[i]=2*pi-angles[i]
    p=zip(angles,points)
    p.sort(key=lambda x:x[0])
    return zip(*p)[1]


#Given a list of points that belong to the same polyhedron and
#the list of hyperplanes that those points are on
#Return a list of polygons that make that polyhedron
def points2polyhedron(points,planes,plotting=True):
    polyhedron=[list() for i in range(len(planes))]
    for i,plane in enumerate(planes):
        cnt=0
        for j,point in enumerate(points):
            a=fabs(dot(point,plane[0:3])+plane[3])
            if a<1e-9:
                cnt+=1
                polyhedron[i].append(point)
        if cnt==0:
            print "ERROR: points2polyhedron, hyperplane with no vertices... something is wrong here."
            exit(0)
        #Need to sort out the polygons if you plan on plotting these
        if plotting:
            polyhedron[i]=points2polygon(polyhedron[i],plane)
    return polyhedron
