
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
    if sum(atoms[:,0])/len(atoms) < 1.0:
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

    print "ERROR: This method (partpaircor) needs to be re-written using the PCPcode shortcut."
    print "Maybe somebody should get off their lazy ass and write it. *cough*"
    exit(0)
    """
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
    """ 
    return [rbins,rdist]

#==================================================================
PCAcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,a,d;
int jn,kn,cn=0;
double dang=180./nbins;

double pi2=2.0*3.14159265;
for(int i=0; i<inloop; i++){
    ajx=atoms[i*3];
    ajy=atoms[i*3+1];
    ajz=atoms[i*3+2];
    for(int j=0;j<nneighbsf[i];j++){
        jn=neighbsf[cn+j];
        aix=atoms[jn*3];
        aiy=atoms[jn*3+1];
        aiz=atoms[jn*3+2];

        //Periodic Bounds
        d = aix-ajx;
        if(d>l[0]/2.0) aix -= l[0];
        if(d<-l[0]/2.0) aix += l[0];
        d = aiy-ajy;
        if(d>l[1]/2.0) aiy -= l[1];
        if(d<-l[1]/2.0) aiy += l[1];
        d = aiz-ajz;
        if(d>l[2]/2.0) aiz -= l[2];
        if(d<-l[2]/2.0) aiz += l[2];

        if(i==jn) 
          continue;
        for(int k=0; k<nneighbsf[i];k++){
            kn=neighbsf[cn+k];
            if(i==kn || kn==jn) 
              continue;
            akx=atoms[kn*3];
            aky=atoms[kn*3+1];
            akz=atoms[kn*3+2];

            //Periodic Bounds
            d = akx-ajx;
            if(d>l[0]/2.0) akx -= l[0];
            if(d<-l[0]/2.0) akx += l[0];
            d = aky-ajy;
            if(d>l[1]/2.0) aky -= l[1];
            if(d<-l[1]/2.0) aky += l[1];
            d = akz-ajz;
            if(d>l[2]/2.0) akz -= l[2];
            if(d<-l[2]/2.0) akz += l[2];

            //Calculate Angle
            dij = (aix-ajx)*(aix-ajx) + (aiy-ajy)*(aiy-ajy) + (aiz-ajz)*(aiz-ajz);
            dik = (aix-akx)*(aix-akx) + (aiy-aky)*(aiy-aky) + (aiz-akz)*(aiz-akz);
            djk = (ajx-akx)*(ajx-akx) + (ajy-aky)*(ajy-aky) + (ajz-akz)*(ajz-akz);  
            x=(dij + djk - dik)/(2.0*sqrt(dij)*sqrt(djk));
            if(fabs(fabs(x)-1.0) <= 1e-9)
              a=0.0;
            else
              a=(360*(acos(x)/pi2+1.0));
            a-= static_cast<double>( static_cast<int>( a / 180.0 ) ) * 180.0;
            bins[(int)(a/dang)]+=1;
        }
    }
    cn+=nneighbsf[i];
}
"""

def paircor_ang(atoms,neighbs,basis,inloop=0,nbins=360,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180

    bins = zeros(nbins)

    if inloop==0:
        inloop=len(atoms)
    
    natoms=len(atoms)
    atoms.shape=natoms*3
    nneighbsf=array([len(i) for i in neighbs])
    neighbsf=array([i for i in flatten(neighbs)])
    l=array([basis[0][0],basis[1][1],basis[2][2]])
    weave.inline(PCAcode,['atoms','neighbsf','nneighbsf','bins','nbins','inloop','l'])
    atoms.shape=[len(atoms)/3,3]
    
    """
    ascale = 180.0/nbins
    for i in range(inloop):
        for ind,j in enumerate(neighbs[i]):
            for k in neighbs[i]:
                if k==j or j==i or k==i: continue
                a = (degrees(ang(atoms[j],atoms[i],atoms[k]))+360.0)%180.0
                bins[int(a/ascale)]+=1
    """
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    
    return [abins,bins]
    
