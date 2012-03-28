#!/usr/bin/python
from math import *

def float2rgb(val, cmin, cmax):
    x = float(val-cmin)/float(cmax-cmin)
    r = min(max(4*fabs(x-0.5)-1., 0.), 1.)
    g = min(max(4*(x-0.25), 0.), 1.)
    b = min(max(4*(0.75-x), 0.), 1.)
    return (r,g,b)

#Blue-Green-Yellow-Red-Magenta
def float2rgb_all(val,cmin,cmax):
    x = (val-cmin)/(cmax-cmin)
    r = min(max( 4*x            -1.0, 0.), 1.)
    g = min(max(-4*fabs(0.375-x)+1.5, 0.), 1.)
    b = min(max( 4*fabs(0.500-x)-1.0, 0.), 1.)
    return (r,g,b)

#white-yello-red-black
def float2rgb_fire(val,cmin,cmax):
    x = (val-cmin)/(cmax-cmin)
    r = min(max(2.0-1.8*x, 0.), 1.)
    g = min(max(1.5-2.0*x, 0.), 1.)
    b = min(max(1.0-2.0*x, 0.), 1.)
    return (r,g,b)
