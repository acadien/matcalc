#!/usr/bin/python

#A collection of helper functions and general utilities for analzing a set of atoms

#from math import *
from scipy import weave,special
from scipy.weave import converters
from numpy import *
import pylab as pl
import numpy as np
#mine
from datatools import flatten,windowAvg 

#=========== Create Ghost Atoms =========================
makeghostcode = """
//Assumes orthogonal basis

double lx = basis[0];
double ly = basis[4];
double lz = basis[8];
double bndlx = -6.0, bndhx = lx+6.0;
double bndly = -6.0, bndhy = ly+6.0;
double bndlz = -6.0, bndhz = lz+6.0;
double ax,ay,az,nax,nay,naz;

int curr = (int)na;
for(int i=0;i<(int)na;i++){

    ax=allatoms[i*3];
    ay=allatoms[i*3+1];
    az=allatoms[i*3+2];

    for(int t1=-1;t1<2;t1++){
        nax = ax+t1*lx;

        if(nax<bndhx && nax>=bndlx)
        for(int t2=-1;t2<2;t2++){
            nay = ay+t2*ly;

            if(nay<bndhy && nay>=bndly)
            for(int t3=-1;t3<2;t3++){
                if(t1==0 && t2==0 && t3==0)
                    continue; //skip self cell

                naz = az+t3*lz;
                if(naz<bndhz && naz>=bndlz){
                    allatoms[curr*3]=nax;
                    allatoms[curr*3+1]=nay;
                    allatoms[curr*3+2]=naz;
                    ghostReduce[curr]=i;
                    curr++;
                }
    }}}

}

return_val=curr;
"""
def makeGhosts(atoms,basis):
    #Creates a layer of ghost atoms around your super cell
    #ghostReduce is the index of the real atom corresponding to the ghost atom

    #Prepare and reshape for weave
    na = len(atoms) #counter of atoms + ghosts
    atoms.shape = [na*3]
    allatoms = zeros([len(atoms)*3*27]) 
    ghostReduce = zeros([len(atoms)*27])
    allatoms[:na*3] = atoms
    ghostReduce[:na] = range(na)
    basis.shape = [9]

    nall = weave.inline(makeghostcode,['na','allatoms','basis','ghostReduce'])

    #reshape for python
    basis.shape = [3,3]
    atoms.shape = [len(atoms)/3,3]
    allatoms=allatoms[:nall*3]
    allatoms.shape = [nall,3]
    ghostReduce = ghostReduce[:nall]

    return allatoms,ghostReduce

#===========  Helper Geometry Functions  ================
#returns the scalar distance between the 2 lists/arrays p1 and p2
#1.742/1.676
#1.65/1.55
distcode = """
double c=0.0;
for(int i=0;i<3;i++)
  c = c + (a[i]-b[i])*(a[i]-b[i]);
return_val=sqrt(c);
"""
def dist(a,b):
    return weave.inline(distcode,['a','b'])

distcodeP = """
double c=0.0,d=0.0;
for(int i=0;i<3;i++){
  d = a[i]-b[i];
  if(d>lengths[i]/2.0) d -= lengths[i];
  if(d<-lengths[i]/2.0) d += lengths[i];
  c = c + d*d;
}
return_val=sqrt(c);
"""
#Calculates dist between 2 points enforcing periodic boundary conditions on orthogonal unit cell
def dist_periodic(a,b,lengths):
    return weave.inline(distcodeP,['a','b','lengths'])


translation_unit_vectors = array(zip( *[ [-1]*9 + [0]*9 + [1]*9       ,\
                                     ( [-1]*3 + [0]*3 + [1]*3 )*3   ,\
                                       [-1,0,1] * 9                ]))

#Calculates dist between 2 points enforcing periodic boundary conditions in 3D, any basis
def minImageDist(a,b,basis):
    #Find all basis vector translations
    translations = tensordot(basis.T, translation_unit_vectors,[1,1]).T
    imageDists = (b-a)+translations

    #Calculate all distances, return the smallest
    return sqrt(min(diag(tensordot(imageDists,imageDists.T,1))))

def minImageAtom(a,b,basis):
    translations = tensordot(basis.T, translation_unit_vectors,[1,1]).T
    imageDists = (b-a)+translations
    ds = diag(tensordot(imageDists,imageDists.T,1))
    bnew = b+ translations[where(ds==min(ds))[0][0]]
    return bnew

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

angdegcode = """
double ang,x,dab=0.0,dac=0.0,dbc=0.0;
double pi2=2.0*3.14159265;

for(int i=0;i<3;i++)
  dab = dab + (a[i]-b[i])*(a[i]-b[i]);
for(int i=0;i<3;i++)
  dac = dac + (a[i]-c[i])*(a[i]-c[i]);
for(int i=0;i<3;i++)
  dbc = dbc + (b[i]-c[i])*(b[i]-c[i]);
x=(dab + dbc - dac)/(2.0*sqrt(dab)*sqrt(dbc));
if(fabs(fabs(x)-1.0) <= 1e-9)
  return_val=0.0;
else {
  ang=360*(acos(x)/pi2+1.0);
  return_val = ang - static_cast<double>( static_cast<int>( ang / 180.0 ) ) * 180.0;
}
"""
def ang_deg(a,b,c):
    return weave.inline(angdegcode,['a','b','c'])

#For spherical angles (theta,phi)
sphang_code="""
double x=b[0]-a[0];
double y=b[1]-a[1];
double z=b[2]-a[2];
double r=sqrt(x*x+y*y+z*z);
angs[0]=acos(z/r); //theta
angs[1]=atan(y/x); //phi
"""
def sphang(a,b):
    angs=zeros(2)
    weave.inline(sphang_code,['a','b','angs'])
    if angs[1] != angs[1]: #test if phi is NaN, occurs if a and b lay on z-axis
        angs[1]=0
    return angs[0],angs[1] #theta,phi

#the spherical harmonic corresponding to two atoms a and b
def pairSphereHarms(a,b,l):
    theta,phi=sphang(a,b)
    return special.sph_harm(np.array(range(-l,l+1)),l,theta,phi)

bounds_code="""
double d=0.0;
for(int i=0;i<3;i++){
  d = a[i]-b[i];
  if(d>lengths[i]/2.0) c[i] = d[i]-lengths[i];
  else if(d<-lengths[i]/2.0) c[i] = d[i]+lengths[i];
  else c[i] = d[i];
}
"""
#Returns the point which is 'b' closest to 'a' employing the bounds given.
def applybounds(a,b,lengths):
    c=zeros(3)
    weave.inline(bounds_code,['a','b','c','lengths'])

def mag(vec):
    return sqrt(dot(vec,vec))

def normalize(vec):
    return vec/mag(vec)

"""
def rotmatx(vec1,vec2):#Returns the rotation matrix that rotates vector1 into vector2
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
"""

def rotmatx(u,v):
    a=normalize(u)
    b=normalize(v)

    w=normalize(cross(a,b))
    U=array([a,w,normalize(cross(a,w))])
    V=array([b,w,normalize(cross(b,w))])
    return dot(V.T,U)

def volume(a,b,c):
    return fabs(dot(a,cross(b,c)))

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
