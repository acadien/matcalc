#!/usr/bin/python

from math import *
from scipy import weave
from scipy.weave import converters

#=========== Efficient implementations of commonly used functions ========

#====================
#  Normal Gaussians
#====================

#f(x)=1/(sig sqrt(2*pi)) e ^ -(x-mu)**2/(2sig^2)
gaussNorm1Dcode = """
double pre=0.3989422804014327; 
double xmus=(x-mu)/sig;
return_val = ( pre/sig ) * exp( xmus*xmus / -2. );
"""

def gaussNorm1D(x,mu,sig):
    return weave.inline(gaussNorm1Dcode,['x','sig','mu'])

#f(x,y)=1/(2pi sigx sigy sqrt(1-rho^2)) e ^ z (z=a mess, check wiki)
gaussNorm2Dcode = """
double xmus=(x-mux)/sigx;
double ymus=(y-muy)/sigy
double rho2=1-rho*rho;
double pre=0.15915494309189535 / (sigx*sigy*sqrt(rho2)); // 1/2pi
return_val = pre * exp( (xmus*xmus + ymus*ymus - 2*rho*xmus*ymus) / (-2.*rho2) );
"""
def gaussNorm2D(x,sigx,mux,y,sigy,muy,rho):
    return weave.inline(gaussNorm2Dcode,['x','y','sigx','sigy','mux','muy','rho'])


#=====================
# Functional Gaussians
#=====================

#f(x)=a e ^ ( (x-x0)^2 / (-2sig^2) ) 
gauss1Dcode = """
double xx=(x-x0)/sig;
return_val = a*exp(xx*xx/-2.);
"""
def gauss1D(a,x,x0,sig):
    return weave.inline(gauss1Dcode,['a','x','x0','sig'])


#f(x,y)=a e ^ ( (x-x0)^2 / ( -2sigx^2) + (y-y0)^2 / ( -2sigy^2) )
gauss2Dcode = """
double xx=(x-x0)/sigx;
double yy=(y-y0)/sigy;
return_val = a * exp( (xx*xx + yy*yy) / -2. );
"""
def gauss2D(a,x,x0,sigx,y,y0,sigy):
    return weave.inline(gauss2Dcode,['a','x','x0','sigx','y','y0','sigy'])
