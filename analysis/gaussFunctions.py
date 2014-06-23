#!/usr/bin/python

from math import *
from scipy import array,zeros
from scipy import weave
from scipy.weave import converters
import numpy

#=========== Efficient implementations of commonly used functions ========

#====================
#  Normal Gaussians
#====================

#f(x)=1/(sig * sqrt(2*pi)) e ^ -(x-mu)**2/(2sig^2)
gaussPointcode = """
double pre=0.3989422804014327; 
double xmus=(x-mu)/sig;
return_val = ( pre/sig ) * exp( xmus*xmus / -2. );
"""

def gaussPoint(x,mu,sig):
    return weave.inline(gaussPointcode,['x','sig','mu'])

#f(x)=1/(sig * sqrt(2*pi)) e ^ -(x-mu)**2/(2sig^2)
gaussNorm1Dcode = """
double pre=0.3989422804014327/sig; 
double xmus;
for(int i=0;i<nx;i++){
    /*
    if(i==0)
        dx = xs[i+1]-xs[i]; 
    if(i==nx-1)
        dx = xs[i]-xs[i-1];
    if(i>0 && i<nx-1)
        dx = (xs[i+1]-xs[i-1])*0.5;
    */
    xmus=(xs[i]-mu)/sig;
    ys[i] = pre * exp( xmus*xmus / -2. )*yv[i]*dx;
}
"""

def gaussNorm1D(xs,mu,sig,yv=None,dx=1.0):
    nx=xs.shape[0]
    ys=zeros(nx)
    if yv==None:
        yv=zeros(nx)+1
    weave.inline(gaussNorm1Dcode,['xs','nx','ys','sig','mu','yv','dx'])
    return ys


#f(x,y)=1/(2pi sigx sigy sqrt(1-rho^2)) e ^ z (z=a mess, check wiki)
gaussNorm2Dcode = """
double xmus=(x-mux)/sigx;
double ymus=(y-muy)/sigy
double rho2=1-rho*rho;
double pre=0.15915494309189535 / (sigx*sigy*sqrt(rho2)); // 1/2pi
return_val = pre * exp( (xmus*xmus + ymus*ymus - 2*rho*xmus*ymus) / (-2.*rho2) );
"""
def gaussNorm2D(xs,ys,params):
    sigx,mux,sigy,muy,rho=params
    for x,y in zip(xs,ys):
        yield weave.inline(gaussNorm2Dcode,['x','y','sigx','sigy','mux','muy','rho'])


#=====================
# Functional Gaussians
#=====================

#f(x)=a e ^ ( (x-x0)^2 / (-2sig^2) ) 
gauss1Dcode = """
double xx;
for(int i=0;i<n;i++){
  xx=( xs[i] - x0 ) / sig;
  ys[i] += a * exp( xx*xx / -2. );
}
"""
def gauss1D(xs,ys,a,x0,sig):
    a,x0,sig=map(float,[a,x0,sig])
    n=len(xs)
    weave.inline(gauss1Dcode,['a','n','x0','sig','xs','ys'])
    return ys

    
#f(x,y)=a e ^ ( (x-x0)^2 / ( -2sigx^2) + (y-y0)^2 / ( -2sigy^2) )
gauss2Dcode = """
double xx=(x-x0)/sigx;
double yy=(y-y0)/sigy;
return_val = a * exp( (xx*xx + yy*yy) / -2. );
"""
def gauss2D(params,xs,ys):
    a,x0,sigx,y0,sigy=map(float,params)
    #print a,x0,sigx,y0,sigy
    return [weave.inline(gauss2Dcode,['a','x','x0','sigx','y','y0','sigy']) for x,y in zip(xs,ys)]
    #return [a*exp(-(((x-x0)/sigx)**2+((y-y0)/sigy)**2)/2.) for x,y in zip(xs,ys)]
    
#Returns the coefs corresponding to the moments of the data
#This is helpful for initializing leastsq fits    
def gauss2Dmoments(xs,ys):
    data=numpy.array(zip(xs,ys))
    total=data.sum()
    X,Y=numpy.indices(data.shape)
    x0 = (X*data).sum()/total
    y0 = (Y*data).sum()/total
    print x0,y0
    col = data[:,int(y0)]
    row = data[int(x0),:]
    sigx=sqrt(abs((numpy.arange(col.size)-y0)**2*col).sum()/col.sum())
    sigy=sqrt(abs((numpy.arange(row.size)-x0)**2*row).sum()/row.sum())
    a=data.max()
    return a,x0,sigx,y0,sigy
    
def gauss1Dmoments(xs,ys):
    x0 = (xs*ys).sum()/ys.sum()
    sigx=sqrt(fabs(sum((xs-x0)**2*ys)/ys.sum()))
    a=ys.max()
    return [a,x0,sigx]
