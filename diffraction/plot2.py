#!/usr/bin/python
import sys
import pylab as pl
import re
import os
import numpy as np

filename=sys.argv[1]
f=open(filename,'r');
i=0;
length=len(f.readlines());
f.close();
f=open(filename,'r');
xs=np.zeros(length);
y1s=np.zeros(length);
y2s=np.zeros(length);
for line in f:
    try:
        [x,y1,y2]=line.split(',')
    except:
        break
    xs[i]=x
    #print x,y1,y2
    y1s[i]=y1
    y2s[i]=y2
    i+=1
pl.subplot(121)
pl.plot(xs[0:i],y1s[0:i],'.')
pl.subplot(122)
pl.plot(xs[0:i],y2s[0:i],'.')
pl.show()
