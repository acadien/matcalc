#!/usr/bin/python

#calculates the dynamic susceptibility of the structure factor
import plotRemote as pr

import sys

#mine
import utils
from interpolate import interp1d
from scipy import fftpack
import numpy as np

utils.usage(["<.isfS file>"],1,1)

isfsFile = sys.argv[1]
if "isfS" not in isfsFile:
    print "Expected .isfS file for input"

logTime = list()
isfs = list()
for isfSline in open(isfsFile,"r").readlines():
    try:
        a = map(float,isfSline.split())
        t = a.pop(0)
        i = sum(a)/len(a)
    except ValueError:
        continue
    logTime.append(t)
    isfs.append(i)

import pylab as pl
pl.semilogx(logTime,isfs)
pr.prshow()
exit(0)
N = len(logTime)
m = max(logTime)
linTime = logTime#[float(i)*m/N for i in range(N)]
linISFs = isfs#interp1d(logTime,isfs,linTime)

T = m/N
xf = np.linspace(0,0.5/T, N/2)
yf = fftpack.fft(linISFs)

pl.semilogx(xf, yf[:N/2])
#pl.semilogx(xf, np.imag(yf[:len(isfs)/2])*xf/len(isfs))
pr.prshow()

