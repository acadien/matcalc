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
def interp3d(ipnt,bnds,points,data):
    interp3dcode="""
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
"""
    [xp,yp,zp]=points
    nps=array(map(len,points))
    return weave.inline(interp3dcode,['ipnt','bnds','data','xp','yp','zp','nps']);

