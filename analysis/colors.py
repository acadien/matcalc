#!/usr/bin/python
from math import *

def float2rgb(val, cmin, cmax):
    x = float(val-cmin)/float(cmax-cmin)
    r = min(max(4*fabs(x-0.5)-1., 0.), 1.)
    g = min(max(4*(x-0.25), 0.), 1.)
    b = min(max(4*(0.75-x), 0.), 1.)
    return (r,g,b)
