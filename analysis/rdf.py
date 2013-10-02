
from math import *
from numpy import degrees,array,zeros
from scipy import weave
from scipy.weave import converters
import pylab as pl
#mine
from struct_tools import volume
from datatools import windowAvg,flatten

#==================================================================
RDFcode="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double dr=cut/nbins;
double cut2=cut*cut;
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=0;j<natoms;j++){
        if(i==j)
          continue;

        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];

        d=aix-ajx;        
        if(d>l[0]/2.0) d -= l[0];
        if(d<-l[0]/2.0) d += l[0];
        c=d*d;

        d = ajy-aiy;
        if(d>l[1]/2.0) d -= l[1];
        if(d<-l[1]/2.0) d += l[1];
        c+=d*d;

        d = ajz-aiz;
        if(d>l[2]/2.0) d -= l[2];
        if(d<-l[2]/2.0) d += l[2];
        c+=d*d;

        if(c<=cut2){
            c=sqrt(c);
            bins[(int)(c/dr)]+=1;
        }
    }
}
"""
def rdfHelper(atoms,basis,bins,cut):
    natoms=len(atoms)
    atoms.shape=len(atoms)*3
    nbins=len(bins)
    l=array([basis[0][0],basis[1][1],basis[2][2]])

    weave.inline(RDFcode,['atoms','natoms','bins','nbins','cut','l'])
    atoms.shape=[len(atoms)/3,3]
    return bins

#==================================================================
def rdf(atoms,basis,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #inloop: number of atoms to do the summation over
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro

    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atoms)
    rdist=rdfHelper(atoms,basis,rdist,cutoff)
    rbins=[(i+0.5)*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    dr=float(cutoff)/nbins
    #Ndensity=N/volume(basis)
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        rdist[i]/=vol

    return [rbins,rdist]

#==================================================================
RDFPerCodeOld="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double dr=cut/nbins;
double cmin;
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=i+1;j<natoms;j++){
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];
        cmin=100000.;
        //Minimum image distance
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
            bins[(int)(cmin/dr)]+=2;
    }
}
"""
RDFPerCode="""
double aix,ajx,aiy,ajy,aiz,ajz,c,d;
double dr=cut/nbins;
double cmin;
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=i+1;j<natoms;j++){
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];
        cmin=100000.;
        //Minimum image distance
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

            if(c<=cut)
                bins[(int)(c/dr)]+=2;
        }}}
    }
}
"""
def rdfperHelper(atoms,bins,cut,b):
    natoms=len(atoms)
    atoms.shape=len(atoms)*3
    b.shape=9
    nbins=len(bins)
    weave.inline(RDFPerCode,['atoms','natoms','bins','nbins','cut','b'])
    atoms.shape=[len(atoms)/3,3]
    b.shape=[3,3]
    return bins


def rdf_periodic(atoms,basis,cutoff=10.0,nbins=1000):
    #atoms: list of atoms[N][3]
    #cutoff: float, max radius to measure radial distro out to
    #nbins: number of bins to store in radial distro
    basis=array(basis)
    atoms=array(atoms)

    bt=basis.T

    if sum(atoms[:,0])/len(atoms) < 1.0:
        atomsp=array([bt.dot(atom) for atom in atoms])
    else:
        atomsp=array(atoms)

    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atomsp)
    rdist=rdfperHelper(atomsp,rdist,cutoff,basis)
    rbins=[i*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    Ndensity=N/volume(basis)
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        rdist[i]/=vol*(Ndensity*N)

    return [rbins,rdist]

#==================================================================

def sf(atoms,basis,cutoff=12.0,nbins=1000):
    basis=array(basis)
    atoms=array(atoms)

#==================================================================
ADFcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,a,d;
int jn,kn,cn=0,cnt=0;
double dang=180./nbins;

double djx,djy,djz,dkx,dky,dkz,djr,dkr,c;

double pi=3.14159265;
for(int i=0; i<natoms; i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=0;j<nneighbsf[i];j++){
        jn=neighbsf[cn+j];
        if(i<=jn) 
          continue;
        
        ajx=atoms[jn*3];
        ajy=atoms[jn*3+1];
        ajz=atoms[jn*3+2];

        //Periodic Bounds
        d = ajx-aix;
        if(d>l[0]/2.0) ajx -= l[0];
        if(d<-l[0]/2.0) ajx += l[0];
        d = ajy-aiy;
        if(d>l[1]/2.0) ajy -= l[1];
        if(d<-l[1]/2.0) ajy += l[1];
        d = ajz-aiz;
        if(d>l[2]/2.0) ajz -= l[2];
        if(d<-l[2]/2.0) ajz += l[2];

        djx = ajx-aix;
        djy = ajy-aiy;
        djz = ajz-aiz;
        djr = djx*djx+djy*djy+djz*djz;

        for(int k=j+1; k<nneighbsf[i];k++){
            kn=neighbsf[cn+k];
            //if(i==kn || kn==jn) 
            //  continue;

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

            dkx = akx-aix;
            dky = aky-aiy;
            dkz = akz-aiz;
            dkr = dkx*dkx+dky*dky+dkz*dkz;

            x = (djx*dkx + djy*dky + djz*dkz)/sqrt(djr*dkr);
            a = acos(x)*180.0/pi;

            //Calculate Angle
            //dij = (aix-ajx)*(aix-ajx) + (aiy-ajy)*(aiy-ajy) + (aiz-ajz)*(aiz-ajz);
            //dik = (aix-akx)*(aix-akx) + (aiy-aky)*(aiy-aky) + (aiz-akz)*(aiz-akz);
            //djk = (ajx-akx)*(ajx-akx) + (ajy-aky)*(ajy-aky) + (ajz-akz)*(ajz-akz);  
            //x=(dij + djk - dik)/(2.0*sqrt(dij)*sqrt(djk));//possibly switch ik - jk
            //if(fabs(fabs(x)-1.0) <= 1e-9)
            //  a=0.0;
            //else
            //  a=180*acos(x)/pi;
            //Bin the bond angle
            //if(a<=180)
            bins[(int)(a/dang)]+=0.5;
            //else
            //  c=dky;
        }
    }
    cn+=nneighbsf[i];
    return_val=c;
}
"""

#angular distribution function
def adf(atoms,neighbs,basis,nbins=360,angtype='deg'):
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
    weave.inline(ADFcode,['atoms','natoms','neighbsf','nneighbsf','bins','nbins','l'],compiler=('gcc'))
    atoms.shape=[len(atoms)/3,3]    
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    bins /= nbins

    return [abins,bins]
    

#==================================================================
RDFBYADFcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,a,d;
double bondij,bondjk;
int jn,kn,cn=0;
double dang=180./nbins;
int binIndex;

double pi2=2.0*3.14159265;
/*
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
*/
}
"""

def rdf_by_adf(atoms,neighbs,basis,nbins=360,angtype='deg'):
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
    weave.inline(RDFBYADFcode,['atoms','natoms','neighbsf','nneighbsf','bins','nbins','MAXNeighbs','bcounts','l'])
    atoms.shape=[len(atoms)/3,3]
    
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    
    return [abins,bins]

#Generates a cutoff based on the RDF of a collection of atoms
def generateRCut(atoms,basis,debug=False):
    #Set rcut to be the first minimum of g(r)
    rvals,gr = rdf_periodic(atoms,basis,cutoff=6.0)

    #Smoothed G(r)
    sgr=windowAvg(gr,n=25).tolist()
    #derivative of smoothed-G(r)
    dsgr = windowAvg(windowAvg([(sgr[i+1]-sgr[i])/(rvals[1]-rvals[0]) for i in range(len(sgr)-1)],n=50),n=20).tolist()#[47:]+[0]*47
    
    #Find the first minima by searching for the first 2 maxima and finding the minima between the two
    #More robust version uses first point at which the first derivative becomes positive after the first peak
    first_neg = [i for i,v in enumerate(dsgr) if v<0][0]
    first_peak = sgr.index(max(sgr[:first_neg]))
    rindex = first_neg + [i for i,v in enumerate(dsgr[first_neg:]) if v>=0][0]
    rcut = rvals[rindex]

     #Should probably check out this plot before continuing
    if debug:
        print rcut
        pl.plot(rvals,[i for i in sgr],lw=3,c="green",label="Smooth G(r)")
        pl.plot(rvals[1:],dsgr,c="red",label="Smooth G'(r)")
        pl.plot([rcut,rcut],[min(gr),max(gr)],lw=3,c="black",label="rcut")
        pl.plot(rvals,gr,c="blue",label="G(r)")
        pl.legend(loc=0)
        pl.show()
        
    return rcut
