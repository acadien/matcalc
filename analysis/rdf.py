
from math import *
from numpy import degrees,array,zeros,linalg,ceil
from scipy import weave
from scipy.weave import converters
import pylab as pl
#mine
from struct_tools import volume,minImageDist
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
for(int i=0;i<natoms;i++){
    aix=atoms[i*3];
    aiy=atoms[i*3+1];
    aiz=atoms[i*3+2];
    for(int j=i+1;j<natoms;j++){
        ajx=atoms[j*3];
        ajy=atoms[j*3+1];
        ajz=atoms[j*3+2];

        //Minimum image distance
        for(int t1=-tmax[0]+1;t1<tmax[0];t1++){
        for(int t2=-tmax[1]+1;t2<tmax[1];t2++){
        for(int t3=-tmax[2]+1;t3<tmax[2];t3++){
            d=aix-ajx+t1*b[0]+t2*b[3]+t3*b[6];
            c=d*d;
            d=aiy-ajy+t1*b[1]+t2*b[4]+t3*b[7];
            c+=d*d;
            d=aiz-ajz+t1*b[2]+t2*b[5]+t3*b[8];
            c+=d*d;

            //if(c==0.0) continue;

            c=sqrt(c);
            if(c<=cut)
                bins[(int)(c/dr)]+=2;
        }}}
    }
}
"""
def rdfperHelper(atoms,bins,cut,b):

    #align things for passing off to weave
    natoms=len(atoms)
    atoms.shape=len(atoms)*3

    #Handles periodic minimum distance atoms to the proper cutoff
    ls=[linalg.norm(bv) for bv in b]
    tmax=array([max([ceil(cut/l*2.0),2]) for l in ls])

    #more aligning
    b.shape=9
    nbins=len(bins)

    weave.inline(RDFPerCode,['atoms','natoms','bins','nbins','cut','b','tmax'])
    
    #re-align back to original shape
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

    rdist=zeros(nbins)
    dr=float(cutoff)/nbins
    N=len(atoms)
    rdist=rdfperHelper(atoms,rdist,cutoff,basis)
    rbins=[i*dr for i in range(nbins)] #the central point of each bin (x-axis on plot)

    Ndensity=N/volume(basis)
    for i,r in enumerate(rbins):
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        rdist[i]/=vol*Ndensity*N

    return [rbins,rdist]

#==================================================================

def sf(atoms,basis,cutoff=12.0,nbins=1000):
    basis=array(basis)
    atoms=array(atoms)

#==================================================================
ADFcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,d;
double a;
int jn,kn,cn=0,cnt=0;

double djx,djy,djz,dkx,dky,dkz,djr,dkr;
double cut2=cut*cut;
double pi=3.14159266; //intentionally overestimate pi by a tiny bit to correct acos behavior
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
        for(int t1=-1;t1<2;t1++){
        for(int t2=-1;t2<2;t2++){
        for(int t3=-1;t3<2;t3++){
            djx = ajx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
            djy = ajy + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
            djz = ajz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

            djr = djx*djx+djy*djy+djz*djz;

            if( djr <= cut2)
               goto breakout;
        }}}
        continue;
        breakout:;

        for(int k=j+1; k<nneighbsf[i];k++){
            kn=neighbsf[cn+k];

            akx=atoms[kn*3];
            aky=atoms[kn*3+1];
            akz=atoms[kn*3+2];

            //Minimum image distance
            for(int t1=-1;t1<2;t1++){
            for(int t2=-1;t2<2;t2++){
            for(int t3=-1;t3<2;t3++){
                dkx = akx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
                dky = aky + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
                dkz = akz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

                dkr = dkx*dkx+dky*dky+dkz*dkz;

                if( dkr > cut2)
                    continue;
                if(djr == 0.0 || dkr == 0.0)
                    continue;

                //Calculate Angle
                x = (djx*dkx + djy*dky + djz*dkz)/sqrt(djr*dkr);
                a = acos(x)/pi;
                if(a != a) //acos returns with NAN for invalid x values, catch these and drop them.
                    continue;
                bins[(int)(a*(nbins-1))]+=0.5;
            }}}
        }
    }
    cn+=nneighbsf[i];
}
return_val=cn;
"""

#angular distribution function
def adf(atoms,neighbs,basis,cutoff,nbins=360,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180

    bins = zeros(nbins)
    
    cut=cutoff
    natoms=len(atoms)
    atoms.shape=natoms*3
    b=basis
    b.shape=9
    nneighbsf=array([len(i) for i in neighbs])
    neighbsf=array([i for i in flatten(neighbs)])
    weave.inline(ADFcode,['atoms','natoms','neighbsf','nneighbsf','bins','nbins','b','cut'],compiler=('gcc'))

    b.shape=[3,3]
    atoms.shape=[len(atoms)/3,3]    
    abins = [(i+0.5)*180./nbins for i in range(nbins)]
    bins /= nbins

    return [abins,bins]
    

#==================================================================
RDFBYADFcode="""
double aix,ajx,akx,aiy,ajy,aky,aiz,ajz,akz;
double dij,dik,djk,x,d;
double adf,rdf;
double djx,djy,djz,dkx,dky,dkz,djr,dkr;

int jn,kn,cn=0;

double cut2=cut*cut;
double pi=3.14159266; //intentionally overestimate pi by a tiny bit to correct acos behavior

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
        for(int t1=-1;t1<2;t1++){
        for(int t2=-1;t2<2;t2++){
        for(int t3=-1;t3<2;t3++){
            djx = ajx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
            djy = ajy + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
            djz = ajz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

            djr = djx*djx+djy*djy+djz*djz;

            if( djr <= cut2)
                goto breakout;
        }}}
        continue;
        breakout:;

        for(int k=j+1; k<nneighbsf[i];k++){
            kn=neighbsf[cn+k];

            akx=atoms[kn*3];
            aky=atoms[kn*3+1];
            akz=atoms[kn*3+2];

            //Minimum image distance
            for(int t1=-1;t1<2;t1++){
            for(int t2=-1;t2<2;t2++){
            for(int t3=-1;t3<2;t3++){
                dkx = akx + t1*b[0]+t2*b[3]+t3*b[6] - aix;
                dky = aky + t1*b[1]+t2*b[4]+t3*b[7] - aiy;
                dkz = akz + t1*b[2]+t2*b[5]+t3*b[8] - aiz;

                dkr = dkx*dkx+dky*dky+dkz*dkz;

                if( dkr > cut2)
                    continue;
                if(djr == 0.0 || dkr == 0.0)
                    continue;

                //Calculate Angle
                x = (djx*dkx + djy*dky + djz*dkz)/sqrt(djr*dkr);
                adf = acos(x)/pi;
                if(adf != adf) //acos returns with NAN for bizarre reasons, catch these and drop them.
                    continue;

                //Bin the bond length by the angle
                rdf=sqrt(dkr)/cut;

                bins[(int)(adf*(nADFbins-1))*nRDFbins + (int)(rdf*(nRDFbins-1))]+=0.5;
            }}}
        }

    }
    
    cn+=nneighbsf[i];
}

return_val=0;
"""

def rdf_by_adf(atoms,neighbs,basis,rcut=10.0,nADFbins=45,nRDFbins=100,angtype='deg'):
    #atoms: list of atoms[N][3]
    #neighbs: the *full* neighbor list for atoms
    #angtype: 'deg' or 'rad'
    #nbins: number of bins to store in radial distro
    #The angular range is always 0-180 deg and angles are taken modulo 180
    
    #This function bins bond length by the angles to which they belong.

    bins = zeros(nADFbins*nRDFbins)

    natoms=len(atoms)
    atoms.shape=natoms*3
    nneighbsf=array([len(i) for i in neighbs])
    neighbsf=array([i for i in flatten(neighbs)])
    cut=rcut
    b=basis
    b.shape=9
    weave.inline(RDFBYADFcode,['atoms','natoms','neighbsf','nneighbsf','bins','nADFbins','nRDFbins','cut','b'])
    atoms.shape=[len(atoms)/3,3]
    
    adfVals = [(i+0.5)*180./nADFbins for i in range(nADFbins)]
    rdfVals = [(i+0.5)*rcut/nRDFbins for i in range(nRDFbins)]
    bins.shape=[nADFbins,nRDFbins]
  
    dr=float(rcut)/nRDFbins
    #Ndensity=N/volume(basis)
    for i,r in enumerate(rdfVals):
        m = sum(bins[:,i])
        if i==0:
            vol=4.0*pi*dr*dr*dr/3.0
        else:
            vol=4.0*pi*r*r*dr
        for j in range(nADFbins):
            if m==0:
                m=1
            bins[j,i]/=vol
        
    return [(adfVals,rdfVals),bins]

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
