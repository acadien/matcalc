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
        t,i = map(float,isfSline.split())
    except ValueError:
        continue
    logTime.append(t)
    isfs.append(i)

N = len(logTime)
m = max(logTime)
linTime = [i*m/N for i in range(N)]
linISFs = interp1d(logTime,isfs,linTime)

T = m/N
xf = np.linspace(0,0.5/T, N/2)
yf = fftpack.fft(isfs)

import pylab as pl
pl.plot(xf, np.imag(yf[:N/2])*xf/N)
pl.show()
