#!/usr/bin/python

#Reads in a poscar, calculates S(q) using Baxter's method
import sys
import poscarIO
from orderParam import radialDistribution
from scipy.signal import cspline1d,cspline1d_eval
from scipy.integrate import cumtrapz,simps
from scipy.optimize import leastsq
from numpy import array,dot,cross,pi,floor
import numpy as np
from datatools import superSmooth,wsmooth
import random
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
    ihc = array([simps([hr[abs(ir-it)]*cr[it] for it,t in enumerate(rads)],x=rads) for ir,r in enumerate(rads)])
    return cr + ihc*density

def residual_crhr(coefs,y,x):
    coefs = array([coefs[i] for i in range(len(coefs)) if i*dx<=Rc] + [0 for i in range(len(coefs)) if i*dx>Rc])
    cr = interp(dx,coefs,x)
    hr = y
    hrc = hrc_calc(dr,cr,hr,x)
    """
    print sum((hr-hrc)**2)
    import pylab as pl
    pl.plot(hrc)
    pl.plot(hr)
    xknots=[i for i in range(len(coefs))]
    pl.scatter(xknots,yknots_original,color="red")
    pl.scatter(xknots,coefs)
    pl.show()
    """
    return hr-hrc


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

def sq_calc(rads,cr,qs):
    dx = rads[1]-rads[0]
    cq = array([ 4*pi*simps(rads*cr*np.sin(q*rads)/q,dx=dx) for q in qs])
    return 1./(1. - density*cq)

def residual_qhr(coefs,y,x):
    Q = interp(dx,coefs,x)
    hr = y
    return x*hr - hrq_calc(dr,Q,hr,x)


basis,atypes,atoms,head,poscar = poscarIO.read(open(sys.argv[1]).readlines())
atoms=array(atoms)
basis=array(basis)
rads,gr = radialDistribution(atoms,basis)
rads = array(rads)

hr = wsmooth(gr)[:-1] - 1#superSmooth(rads,gr)[:-20] - 1
#hr = gr - 1
dr=rads[1]-rads[0]
density = atoms.shape[0]/dot(cross(basis[0],basis[1]),basis[2])

rcut = rads[-1]
radsCut = [i for i in rads if i <= rcut]
hrCut = [hr[i] for i in range(len(rads)) if rads[i] <= rcut]

Nknots=70
xknots=array([i*max(radsCut)/Nknots for i in range(Nknots)])
dx=xknots[1]-xknots[0]
yknots=list()
for i in range(Nknots):
    if i==0:
        yknots.append(0)
    else:
        yknots.append(yknots[-1]*0.9+0.1*random.random()-0.05)
yknots=array(yknots)
#yknots=array([-0.5+random.random() for i in range(Nknots)])
#yknots=array([i for i in range(Nknots)])

Rc = 5.0

#(lsq,success) = leastsq(residual_qhr,yknots,args=(hr,rads),maxfev=100)
#Q  = interp(dx,lsq,rads)
#cr = crq_calc(dr,Q,rads)
#hr_new = hrq_calc(dr,Q,hr,rads)
yknots_original = array(list(yknots))
(lsq,success) = leastsq(residual_crhr,yknots,args=(hr,rads),maxfev=1)
cr = interp(dx,lsq,rads)
hr_new = hrc_calc(dr,cr,hr,rads)

Nq = 1000
qMax = 12.0
qs = array([q/float(Nq)*qMax for q in range(Nq)])
sq = sq_calc(rads,cr,qs)

import pylab as pl
pl.plot(rads,hr,label="hrold")
pl.plot(rads,hr_new,label="hrnew")
pl.plot(rads,cr,label="cr")
pl.plot(qs,sq)
#pl.plot(rads,deriv(cr,dr))
#pl.plot([0,rads[-1]],[0,0])
pl.show()
