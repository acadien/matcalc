#!/usr/bin/python
from scipy import *
from scipy import weave
from scipy.weave import converters
from scipy.interpolate import InterpolatedUnivariateSpline

def interp1d(x,y,xi):
    ius = InterpolatedUnivariateSpline(x,y)
    return ius(xi)

#interpolation in 2 variables: x,y
#zxy is how z12 should be read.
def bilinear_interpolation(x, y,(xmin,xmax),(ymin,ymax),(z11,z12,z21,z22)):
    #print z11,z12,z21,z22
    return (z11 * (xmax - x) * (ymax - y) +  z21 * (x - xmin) * (ymax - y) +  \
            z12 * (xmax - x) * (y - ymin) +  z22 * (x - xmin) * (y - ymin) )  \
            / ((xmax - xmin) * (ymax - ymin) )

#3d data interpolate onto points
#interp3d_simp(interpolation pnt, flattened 3D-data, shape of unflattened grid):
def interp3d(ipnt,data,shape):
    interp3dcode="""
int xstd=shape[1]*shape[2]; //x-stride
int ystd=shape[2];     //y-stride

//Cube surrounding requested point
int x1=ipnt[0];
int x2=x1+1;
if( x1 >= shape[0]-1){
    x2 = x1;    
    x1 -= 1;
} 
int y1=ipnt[1];
int y2=y1+1;
if( y1 >= shape[1]-1){
    y2 = y1;    
    y1 -= 1;
}
int z1=ipnt[2];
int z2=z1+1;
if( z1 >= shape[2]-1){
    z2 = z1;    
    z1 -= 1;
}
//First zmin plane
double q11=data[x1*xstd+y1*ystd+z1];
double q12=data[x1*xstd+y2*ystd+z1];
double q21=data[x2*xstd+y1*ystd+z1];
double q22=data[x2*xstd+y2*ystd+z1];
//Bilinear interpolation
double z1val= (q11 * (x2 - ipnt[0]) * (y2 - ipnt[1]) +  q21 * (ipnt[0] - x1) * (y2 - ipnt[1]) +  \
                 q12 * (x2 - ipnt[0]) * (ipnt[1] - y1) +  q22 * (ipnt[0] - x1) * (ipnt[1] - y1) )  \
                / ((x2 - x1) * (y2 - y1));

//Second z2 plane
q11=data[x1*xstd+y1*ystd+z2];
q12=data[x1*xstd+y2*ystd+z2];
q21=data[x2*xstd+y1*ystd+z2];
q22=data[x2*xstd+y2*ystd+z2];
//Bilinear interpolation
double z2val= (q11 * (x2 - ipnt[0]) * (y2 - ipnt[1]) +  q21 * (ipnt[0] - x1) * (y2 - ipnt[1]) +  \
                 q12 * (x2 - ipnt[0]) * (ipnt[1] - y1) +  q22 * (ipnt[0] - x1) * (ipnt[1] - y1) )  \
                / ((x2 - x1) * (y2 - y1));

//return_val=(z1val*fabs(z2-ipnt[2])+z2val*fabs(ipnt[2]-z1))/fabs(z2-z1);
return_val = data[x1*xstd+y1*ystd+z1];
"""
    for i,mx in zip(ipnt,shape):
        if i >= mx or i < 0:
            print "interpolate.interp3d error: Requested point out of bounds %f %f %f"%tuple(ipnt)
            exit(0)

    ipnt=array(ipnt)#for some reason ipnt is completely broken if its not copied
    return weave.inline(interp3dcode,['ipnt','data','shape'])

"""
#3d data interpolate onto points
def interp3d(ipnt,bnds,points,data):
    interp3dcode=
int xstd=nps[1]*nps[2]; //x-stride
int ystd=nps[2];     //y-stride

//Indeces
int xmin=(int)((ipnt[0]-bnds[0])/(bnds[1]-bnds[0])*nps[0]);
int xmax=xmin+1;
int ymin=(int)((ipnt[1]-bnds[2])/(bnds[3]-bnds[2])*nps[1]);
int ymax=ymin+1;
int zmin=(int)((ipnt[2]-bnds[4])/(bnds[5]-bnds[4])*nps[2]);
int zmax=zmin+1;

//Corresponding points
double x1=xp[xmin];
double x2=xp[xmax];
double y1=yp[ymin];
double y2=yp[ymax];
double z1=zp[zmin];
double z2=zp[zmax];

//First zmin plane
double q11=data[xmin*xstd+ymin*ystd+zmin];
double q12=data[xmin*xstd+ymax*ystd+zmin];
double q21=data[xmax*xstd+ymin*ystd+zmin];
double q22=data[xmax*xstd+ymax*ystd+zmin];
//Bilinear interpolation
double z1val= (q11 * (x2 - ipnt[0]) * (y2 - ipnt[1]) +  q21 * (ipnt[0] - x1) * (y2 - ipnt[1]) +  \
                 q12 * (x2 - ipnt[0]) * (ipnt[1] - y1) +  q22 * (ipnt[0] - x1) * (ipnt[1] - y1) )  \
                / ((x2 - x1) * (y2 - y1));

//Second zmax plane
q11=data[xmin*xstd+ymin*ystd+zmax];
q12=data[xmin*xstd+ymax*ystd+zmax];
q21=data[xmax*xstd+ymin*ystd+zmax];
q22=data[xmax*xstd+ymax*ystd+zmax];
//Bilinear interpolation
double z2val= (q11 * (x2 - ipnt[0]) * (y2 - ipnt[1]) +  q21 * (ipnt[0] - x1) * (y2 - ipnt[1]) +  \
                 q12 * (x2 - ipnt[0]) * (ipnt[1] - y1) +  q22 * (ipnt[0] - x1) * (ipnt[1] - y1) )  \
                / ((x2 - x1) * (y2 - y1));

return_val=(z1val*fabs(z2-ipnt[2])+z2val*fabs(ipnt[2]-z1))/fabs(z2-z1);
    [xp,yp,zp]=points
    nps=array(map(len,points))
    return weave.inline(interp3dcode,['ipnt','bnds','data','xp','yp','zp','nps']);

"""
