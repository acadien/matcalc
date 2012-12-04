
from math import *
from numpy import degrees
from scipy import weave
from scipy.weave import converters
import operator
#mine
from struct_tools import *

#==================================================================
PCcode="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double dr=cut/nbins;
for(int i=0;i<inloop;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=0;j<natoms;j++){
        if(i==j) continue;
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];
        d=aix-ajx;
        c=d*d;
        d=aiy-ajy;
        c+=d*d;
        d=aiz-ajz;
        c+=d*d;
        c=sqrt(c);
        if(c<=cut)
            bins[(int)(c/dr)]++;
    }
}
"""
def pairCorHelper(atoms,bins,cut,inloop):
    natoms=len(atoms)
    atoms.shape=len(atoms)*3
    nbins=len(bins)
    weave.inline(PCcode,['atoms','natoms','bins','nbins','cut','inloop'])
    atoms.shape=[len(atoms)/3,3]
    return bins

#==================================================================
def paircor(atoms,inloop=0,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #inloop: number of atoms to do the summation over
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro

    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atoms)
    if inloop==0:
        inloop=N

    rdist=pairCorHelper(atoms,rdist,cutoff,inloop)
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
double cmin;
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=i+1;j<natoms;j++){
        //if(i==j) continue;
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];
        cmin=100000.;
        for(int t1=-1;t1<2;t1++){
        for(int t2=-1;t2<2;t2++){
        for(int t3=-1;t3<2;t3++){
            d=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            c=d*d;
            d=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            c+=d*d;
            d=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            c+=d*d;
            c=sqrt(c);


            if(c<=cmin)
                cmin=c;
        }}}
        if(cmin<cut)
            bins[(int)(cmin/dr)]++;
    }
}
"""
def pairCorPerHelper(atoms,bins,cut,b):
    natoms=len(atoms)
    atoms.shape=len(atoms)*3
    b.shape=9
    nbins=len(bins)
    weave.inline(PCPcode,['atoms','natoms','bins','nbins','cut','b'])
    atoms.shape=[len(atoms)/3,3]
    b.shape=[3,3]
#    bins/=atoms.shape[0]
    return bins


def paircor_periodic(atoms,basis,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro
    basis=array(basis)

    bt=basis.T
    atomsp=atoms

    if sum(atoms[:,0])/len(atoms) < 1.0:
        atomsp=[bt.dot(atom) for atom in atoms]
        #atoms=array(map(lambda x: [x[i]* for i in range(3)],atoms))
    atomsp=array(atomsp)

    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atomsp)
    rdist=pairCorPerHelper(atomsp,rdist,cutoff,basis)
    rbins=[i*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    Ndensity=N/volume(*basis)
    print Ndensity,N
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        rdist[i]/=vol*Ndensity*N/3

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
for(int i=0; i<natoms; i++){
    ajx=atoms[i*3];
    ajy=atoms[i*3+1];
    ajz=atoms[i*3+2];
    for(int j=0;j<nneighbsf[i];j++){
        jn=neighbsf[cn+j];
        if(i==jn) 
          continue;

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

            //Bin the bond angle
            bins[(int)(a/dang)]++;
        }
    }
    cn+=nneighbsf[i];
}
"""

def paircor_ang(atoms,neighbs,basis,nbins=360,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180

    bins = zeros(nbins)

    natoms=len(atoms)
    atoms.shape=natoms*3

    nneighbsf=array([len(i) for i in neighbs])
    neighbsf=array([i for i in flatten(neighbs)])

    l=array([basis[0][0],basis[1][1],basis[2][2]])

    weave.inline(PCAcode,['atoms','natoms','neighbsf','nneighbsf','bins','nbins','l'])

    atoms.shape=[len(atoms)/3,3]    
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    
    return [abins,bins]
    

#==================================================================
PCbyAcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,a,d;
double bondij,bondjk;
int jn,kn,cn=0;
double dang=180./nbins;
int binIndex;

double pi2=2.0*3.14159265;
for(int i=0; i<natoms; i++){
    ajx=atoms[i*3];
    ajy=atoms[i*3+1];
    ajz=atoms[i*3+2];
    for(int j=0;j<nneighbsf[i];j++){
        jn=neighbsf[cn+j];
        if(i==jn) 
          continue;

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

        //Bond length i-j
        bondij=sqrt((aix-ajx)*(aix-ajx)+(aiy-ajy)*(aiy-ajy)+(aiz-ajz)*(aiz-ajz));

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

            //Bond length j-k
            bondjk=sqrt((ajx-akx)*(ajx-akx)+(ajy-aky)*(ajy-aky)+(ajz-akz)*(ajz-akz));

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

            //Bin the bond length by the angle
            binIndex=(int)(a/dang);
            bins[binIndex*MAXNeighbs + bcounts[binIndex]++]=bondij;
            bins[binIndex*MAXNeighbs + bcounts[binIndex]++]=bondjk;
        }
    }
    cn+=nneighbsf[i];
}
"""

def paircor_binByAng(atoms,neighbs,basis,nbins=360,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180
    
    #This function bins bond length by the angles to which they belong.

    MAXNeighbs=100
    bins = zeros(nbins*MAXNeighbs)
    bcounts = zeros(nbins) #extra memory needed for C operation

    natoms=len(atoms)
    atoms.shape=natoms*3
    nneighbsf=array([len(i) for i in neighbs])
    neighbsf=array([i for i in flatten(neighbs)])
    l=array([basis[0][0],basis[1][1],basis[2][2]])
    weave.inline(PCbyAcode,['atoms','natoms','neighbsf','nneighbsf','bins','nbins','MAXNeighbs','bcounts','l'])
    atoms.shape=[len(atoms)/3,3]
    
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    
    return [abins,bins]
