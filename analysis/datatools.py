#!/usr/bin/python

from numpy import *
#from scipy.signal import freqz
from numpy.fft import *
from numpy import convolve
import pylab as pl
from itertools import chain
from scipy import ndimage
from scipy import weave
from scipy.weave import converters
#mine
from gaussFunctions import *

#uses a guassian smooth convoluted with finite differences to get an absurdly smooth line but with edge effects
superSmoothCode="""
double pre=0.3989422804014327/sigma;
double dx,xmus;

for(int a=0;a<N;a++){
    for(int b=0;b<N;b++){
        if(b==0)
            dx = xs[b+1]-xs[b];
        if(b==N-1)
            dx = xs[b]-xs[b-1];
        if(b>1 && b<N-1)
            dx = (xs[b+1]-xs[b-1])/2.0;

        xmus = (xs[a]-xs[b])/sigma;
        smoothys[a] += pre * exp( xmus * xmus * -0.5) * ys[b] * dx;
}}
"""
def superSmooth(xs,ys,sigma=1.0):
    N=len(ys)
    smoothys=zeros(N)
    xs=array(xs)
    ys=array(ys)
    weave.inline(superSmoothCode,['xs','ys','N','smoothys','sigma'])
    return smoothys

#1D data
def gaussSmooth(a,n=9,mode='same'):
    if n%2==0:
        n+=1
    center=(n-1)/2.
    gWeight=[gaussNorm1D(float(i),center,float(1)) for i in range(n)]
    return convolve(a,gWeight,mode=mode)#[center-1:-center-1]

#1D data
def windowAvg(a,n=11):
    #a: the list/array to run the window average over
    #n: the size of the window
    return convolve(a, ones(n)/n,'same')

#This function nabbed from www.scipy.org/Cookbook/SignalSmooth
#Applies spectral methods (various convolution functions) to a window smoothing
def wsmooth(x,window_len=11,window='hanning'):
    if type(x) != type(array(1)):
        x=array(x)
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')

    y=convolve(w/w.sum(),s,mode='valid')
    M=window_len/2
    return y[M:-M+1]

"""
#Finite Impulse Response
def FIR(n,wl,wh):
    return (sin(wh*n) - sin(wl*n))/pi/n

def bandpass(x,y,bandlow,bandhigh):
    #x/y: the lists/arrays to perform the bandpass on, assumes evenly spaced
    #bandlow/bandhigh: the frequency range to work in

    N=len(y)
    dx=x[1]-x[0]

    y=[sin(i) for i in x]
    wl=2*dx*bandlow*pi
    wh=2*dx*bandhigh*pi

    M=N
    bn = [(wh-wl)/pi]
    bn += [FIR(i-M/2,wl,wh) for i in range(1,M)]
    bn = array(bn)
    bn = bn*kaiser(M,5.2)
    

    [w,h]=freqz(bn,1)

    filtered=convolve(bn,array(y))
    pl.plot(y)
    pl.plot(filtered)
    pl.show()
    exit()
"""

#Takes a 2D list and turns it into a 1D list
def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)
