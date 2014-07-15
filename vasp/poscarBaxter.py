#!/usr/bin/python

#Reads in a poscar, calculates S(q) using Baxter's method
import sys
import poscarIO
from orderParam import radialDistribution
from scipy.signal import cspline1d,cspline1d_eval
from scipy.integrate import cumtrapz,simps
from scipy.optimize import leastsq
from numpy import array,dot,cross,pi,floor
from datatools import superSmooth,wsmooth

import pylab as pl

def deriv(xs,h):
    N=len(xs)
    return array([(xs[1]-xs[0])/h] + [(xs[i+1]-xs[i-1])*0.5/h for i in range(1,N-1)] + [(xs[N-1]-xs[N-2])/h])
    
# dx=space between knots, xis=points to interpolate on
def interp(dx,knots,xis):
    cspl = cspline1d(knots) #cupic spline coefficients
    return cspline1d_eval(cspl,xis,dx=dx,x0=0.0)

#calculate h(r) from c(r)
def hrc_calc(dr,cr,hr,rads):
    ihc = array([simps([hr[abs(floor((r-t)/dr))]*cr[it] for it,t in enumerate(rads)],x=rads) for ir,r in enumerate(rads)])
    return cr+ihc*density

#calculate h(r) from Q(r)
def hrq_calc(dr,Q,hr,rads):
    N=len(Q)

    #first derivative from finite differences
    dQ = deriv(Q,dr)

    #integral of Q convoluted with h and r-t
    iQ = array([simps([(r-t)*hr[abs(floor((r-t)/dr))]*Q[it] for it,t in enumerate(rads)],x=rads) for ir,r in enumerate(rads)])

    return (-dQ + 2.0*pi*density*iQ)

#calculate c(r) from Q
def crq_calc(dr,Q,rads):
    N=len(Q)

    #first derivative from finite differences
    dQ = deriv(Q,dr)

    #integral of Q convoluted with its derivative
    iQ = array([simps([Q[abs(ir-it)]*dQ[it] for it,t in enumerate(rads[ir:])],x=rads[ir:]) for ir,r in enumerate(rads)])

    return (-dQ + 2.0*pi*density*iQ)/rads


def sq_calc(cr):
    return

def residual_qhr(coefs,y,x):
    Q = interp(dx,coefs,x)
    hr = y
    return x*hr - hrq_calc(dr,Q,hr,x)

def residual_crhr(coefs,y,x):
    cr = interp(dx,coefs,x)
    hr = y
    return hr-hrc_calc(dr,cr,hr,x)

basis,atypes,atoms,head,poscar = poscarIO.read(open(sys.argv[1]).readlines())
atoms=array(atoms)
basis=array(basis)
rads,gr = radialDistribution(atoms,basis)
hr = wsmooth(gr) - 1#superSmooth(rads,gr)[:-20] - 1
#hr = gr - 1
dr=rads[1]-rads[0]
density = dot(cross(basis[0],basis[1]),basis[2])/atoms.shape[0]

rcut = rads[-1]
radsCut = [i for i in rads if i <= rcut]
hrCut = [hr[i] for i in range(len(rads)) if rads[i] <= rcut]

Nknots=20
xknots=array([i*max(radsCut)/Nknots for i in range(Nknots)])
dx=xknots[1]-xknots[0]
yknots=array([i for i in range(Nknots)])

#(lsq,success) = leastsq(residual_qhr,yknots,args=(hr,rads),maxfev=100)
#Q  = interp(dx,lsq,rads)
#cr = crq_calc(dr,Q,rads)
#hr_new = hrq_calc(dr,Q,hr,rads)

(lsq,success) = leastsq(residual_crhr,yknots,args=(hr,rads),maxfev=100)
cr = interp(dx,lsq,rads)
hr_new = hrc_calc(dr,cr,hr,rads)

import pylab as pl
pl.plot(rads,hr)
pl.plot(rads,hr_new)
pl.plot(rads,cr)
pl.plot(rads,deriv(cr,dr))
pl.plot([0,rads[-1]],[0,0])
pl.show()
