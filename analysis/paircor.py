
from math import *
from numpy import degrees
from scipy import weave
from scipy.weave import converters
import operator
#mine
from struct_tools import *

#==================================================================
def paircor(atoms,inloop=0,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #inloop: number of atoms to do the summation over
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro

    rdist=[0]*nbins
    dr=float(cutoff)/nbins
    N=len(atoms)
    if inloop==0:
        inloop=N

    for i in range(inloop):
        ai=atoms[i]
        for j in range(i+1,N):
            aj=atoms[j]
            r=dist(ai,aj)
            if r>cutoff:
                continue
            rdist[int(r/dr)]+=1
    rbins=[(i+0.5)*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    for i,rad in enumerate(rbins):
        if i==0:
            vol=4.0*pi*rad*rad*rad/3.0
        else:
            vol=4.0*pi*rad*rad*dr
        rdist[i]/=vol
        
    return [rbins,rdist]

#==================================================================
PCPcode="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double dr=cut/nbins;
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=0;j<natoms;j++){
        if(i==j) continue;
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];
        for(double tx=-1.0;tx<2.0;tx+=1.){
        for(double ty=-1.0;ty<2.0;ty+=1.){
        for(double tz=-1.0;tz<2.0;tz+=1.){
            d=aix-ajx+tx*l[0];
            c=d*d;
            d=aiy-ajy+ty*l[1];
            c+=d*d;
            d=aiz-ajz+tz*l[2];
            c+=d*d;
            c=sqrt(c);

            if(c<=cut)
                bins[(int)(c/dr)]++;
        } } }
    }
}
"""
def paircorHelper(atoms,bins,cut,l):
    natoms=len(atoms)
    atoms.shape=len(atoms)*3
    nbins=len(bins)
    weave.inline(PCPcode,['atoms','natoms','bins','nbins','cut','l'])
    atoms.shape=[len(atoms)/3,3]
    return bins

def paircor_periodic(atoms,lengths,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro

    atoms=array(map(lambda x: [x[i]*lengths[i] for i in range(3)],atoms))
    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atoms)

    rdist=paircorHelper(array(atoms),rdist,cutoff,lengths)
    rbins=[i*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    Ndensity=N/reduce(operator.mul,lengths)
    for i in range(nbins):
        r=float(i)*dr
        if i==0:
            vol=3.0*pi*dr*dr*dr/4.0
        else:
            vol=4.0*pi*r*r*dr
        rdist[i]/=vol*Ndensity*N

    return [rbins,rdist]

#==================================================================
def partpaircor(atoms,types,type1,type2,inloop=0,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #types: list of types[N]
    #type1, type2: which types to compare against
    #inloop: number of atoms to do the summation over
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro

    rdist=[0]*nbins
    dr=float(cutoff)/nbins
    N=len(atoms)
    if inloop==0:
        inloop=N

    for i in range(inloop):
        [xi,yi,zi]=atoms[i]
        if types[i]!=type1: continue
        
        for j in range(i+1,N):
            if types[j]!=type2: continue

            [xj,yj,zj]=atoms[j]
            r=((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)**0.5
            if r>cutoff:
                continue
            rdist[int(r/dr)]+=1
    rbins=[(i+0.5)*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    for i,rad in enumerate(rbins):
        if i==0:
            vol=4.0*pi*rad*rad*rad/3.0
        else:
            vol=4.0*pi*rad*rad*dr
        rdist[i]/=vol
        
    return [rbins,rdist]

#==================================================================
def paircor_ang(atoms,neighbs,inloop=0,nbins=1000,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180

    ascale = 180.0/nbins
    aang = [0]*nbins

    if inloop==0:
        inloop=len(atoms)
    for i in range(inloop):
        for ind,j in enumerate(neighbs[i]):
            for k in neighbs[i][ind+1:]:
                a = (degrees(ang(atoms[i],atoms[j],atoms[k]))+360.0)%180.0
                if a<0.0 or a>180.0:
                    print "Warning, angle out of range."
                aang[int(a/ascale)]+=1

    abins = [(i+0.5)*ascale for i in range(nbins)]

    return [abins,aang]
    
