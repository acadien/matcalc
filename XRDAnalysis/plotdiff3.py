#!/usr/bin/python
import sys
import pylab as pl
import os
import numpy as np

n=len(sys.argv)
fig = pl.figure()
while n>0:
    filename=sys.argv[n-1]
    n-=1
    f=open(filename,'r');
    i=0;
    length=len(f.readlines());
    f.close();
    f=open(filename,'r');
    xs=np.zeros(length);
    y1s=np.zeros(length);
    y2s=np.zeros(length);
    y3s=np.zeros(length);
    for line in f:
        try:
            [x,y1,y2,y3]=line.split(',')
        except:
            break
        xs[i]=x
        y1s[i]=y1
        y2s[i]=y2
        y3s[i]=y3
        i+=1

    ax1=fig.add_subplot(131)
    ax1.plot(xs[0:i],y1s[0:i],'.')

    ax2=fig.add_subplot(132)
    ax2.plot(xs[0:i],y2s[0:i],'.')

    ax3=fig.add_subplot(133)
    ax3.plot(xs[0:i],y3s[0:i],'.')

pl.show()
