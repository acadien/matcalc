#!/usr/bin/python
from numpy import *
from scipy.signal import cspline1d,cspline1d_eval
from scipy.optimize import leastsq
from itertools import chain

#My first real code written in python, so I'm sorry its so C-like :(

ncoef=5

def flatten(listoflists):
    b=list()
    for a in listoflists:
        b+=a
    return b

def localMax(seq,start,delta):
    candidate = False
    m=-10000.0
    base=seq[start]
    found=False
    for i,elem in enumerate(seq[start:]):  
        if elem-base >= delta: found=True
        if elem < base: base=elem
        if elem >= m:
            if not(candidate): base=elem
            m = elem        
            candidate = True
        else:      
            if found: 
                return i-1+start
            m = elem        
            candidate = False 

    return -1    

def localMin(seq,start,delta):
    i = start
    candidate = False
    m=10000.0
    base=seq[start]
    found=False
    for elem in seq[start:]:
        if base-elem >= delta: found=True
        if elem > base: base=elem
        if elem <= m:
            if not(candidate): base=elem
            m=elem
            candidate=True
        else:
            if found:
                return i-1
            m=elem
            candidate=False
        i+=1
    return -1

###################
#looks for a steep drop to something near 0 
#shortly after start of sequence, returns this drop point
def finddrop(seq,start): 
    half=(seq[start]+seq[-1])/2
    (a,b)=(0,0)
    for i in range(start,len(seq)):
        if seq[i]<half:
            a=i
            break
    half=(seq[start]+seq[a])/2
    for i in range(start,len(seq)):
        if seq[i]<half:
            b=i
            break
    return b-(a-b)

def loadMgO(filename):
    mgof=open(filename,'r')
    for line in mgof:
        if line[0]=='#':
            continue
        break
    #First non comment line is temperatures
    line=line.rstrip()
    Ts=[float(i) for i in line.split(' ')]
    ARs=list() #a ratio (MgO FCC lattice constant ratio)
    alen=i=0
    PPs=zeros(0)
    for line in mgof:
        if line[0]=='#':
            continue
        elif alen==0:
            alen=line
            PPs=zeros((int(alen),len(Ts)))
            continue
        vals=line.split(' ')
        ARs.append(float(vals[0])**(1.0/3.0))
        PPs[i]=vals[1:]
        i+=1
    return (ARs,Ts,PPs)

def interp_T(T,Ts,PPs):
    xvals=zeros(1)
    xvals[0]=T
    Ps=list()
    for i in PPs:
        cj=cspline1d(i)
        Ps.append(cspline1d_eval(cj,xvals,dx=500.0,x0=Ts[0]))
    return Ps

def splinefit(xdata,ydata):
    xs=arange(min(xdata),max(xdata),0.1/float(len(xdata)))
    cj=cspline1d(array(ydata))
    return (xs,cspline1d_eval(cj,xs,dx=xdata[1]-xdata[0],x0=xdata[0]))
    
def gaussval(x,coef):
    return (1/(sqrt(2*pi*coef[0])))*exp(-(x-coef[1])**2/(2*(coef[0])))

def gval(x,coef,c):
    return coef[0]*exp(-(x-c)**2/(2*(coef[1])))

def lrntzval(x,coef):
    return coef[1]/(1+(coef[1]*pi*(x-coef[0]))**2)

def glval(x,coef):
    tot=fabs(coef[1])*(coef[0]*gaussval(x,coef[2:4])+(1-coef[0])*lrntzval(x,coef[3:5]))+coef[5]*x+coef[6]
    if coef[0]<0 or coef[0]>1:
        tot=fabs(tot)*-1
    return tot

def glresiduals(coef,y,x):
    return y-glval(x,coef)

def gnval(x,coef,n,consts):
    return sum([gval(x,coef[i*2:(i+1)*2],consts[i]) for i in range(n)])
#return gval(x,coef[:3])+gval(x,coef[3:6])+gval(x,coef[6:])

def gnresiduals(coef,y,x,n,consts): #residual from 3 gaussians
    return y-gnval(x,coef,n,consts)

def glinitguess(xdata,ydata):
    if len(xdata)>15:
        cen=len(xdata)/2    
        start=localMin(ydata[0:cen],0,3)
        if start==-1:
            ylist=ydata[0:cen]
            start=ylist.index(min(ylist))+1
            
        end=localMin(ydata[cen:],0,3)
        if end==-1:
            ylist=ydata[cen:]
            end=cen+ylist.index(min(ylist))+1
        else:
            end+=cen

        d=end-start
        if d<7:
            dl=int(7-d)/2+1
            if start-dl>=0:
                start-=dl
            else:
                start=0
            if end+dl<len(xdata):
                end+=dl
            else:
                end=-1
    else:
        start=0
        end=-1

    (xdata,ydata)=splinefit(xdata[start:end],ydata[start:end])
    cen=argwhere(ydata==max(ydata))[0][0]
    scal=(max(ydata)-min(ydata))
    glrat=0.5
    if max(ydata)/min(ydata)>3:
        glrat=0.2
    if scal<500:
        glrat=0.7

    #Initial guess of coefs
    coefs=[glrat,scal**0.8,0.1,xdata[cen],(max(ydata))**(0.5),(ydata[-1]-ydata[0])/(xdata[-1]-xdata[0]),0]

    #make the linear shift (coef[6]) so that the minimum point is at min(ydata):
    ym=min(ydata)
    ymi=where(ydata==ym)[0][0]
    coefs[6]=ym-glval(xdata[ymi],coefs)
    
    return [xdata,ydata,coefs]

def fiterr(data1,data2):
    return sum([fabs(a-b) for a,b in zip(data1,data2)])

def glfit(xdata,ydata,coefs):
    (lsq,success)=leastsq(glresiduals,coefs,args=(ydata,xdata),maxfev=40000)
    return [glval(xdata,lsq),lsq]

def mglval(x,coef,N):
    runtot=zeros(len(x))
    for i in range(N):
        z=i*ncoef
        runtot+=fabs(coef[z+1])*(coef[z]*gaussval(x,coef[z+2:z+4])+(1-coef[z])*lrntzval(x,coef[z+3:z+5]))
        if coef[z]<0 or coef[z]>1:
            return fabs(runtot)*-1

    runtot+=coef[-2]*x+coef[-1]
    return runtot

def mglresiduals(coef,y,x,N):
    return y-mglval(x,coef,N)

def mglinitguess(xdata,ydata,initcoefs):
    #Setup coefficients
    N=len(initcoefs)
    coefs=list()
    pkht=list()
    for cs in initcoefs:
        pkht.append(cs[4])
        for c in cs:
            coefs.append(c)
        coefs.pop()
        coefs.pop()

    #Fix scaling differences in the peaks to make the fit easier
    pkdif=0
    for i in range(4,len(coefs),ncoef):
        pkdif+=(coefs[i]+coefs[i-ncoef])**2
    pkdif/=N
    if pkdif>1000:
        mxpkind=coefs.index(max([coefs[i] for i in range(4,len(coefs),ncoef)]))
        coefs[mxpkind]=(coefs[mxpkind])**(0.7)
   
    scdif=0
    for i in range(1,len(coefs),ncoef):
        scdif+=(coefs[i]+coefs[i-ncoef])**2
    scdif/=N
    if scdif>20000:
        mxscind=coefs.index(max([coefs[i] for i in range(1,len(coefs),ncoef)]))
        coefs[mxscind]=(coefs[mxscind])**(0.8)
        
    coefs.append((ydata[-1]-ydata[0])/(xdata[-1]-xdata[0]))
    coefs.append(min(ydata))

    (xdata,ydata)=splinefit(xdata,ydata)
    return [xdata,ydata,coefs]

def mglfit(xdata,ydata,coefs,N):
    (lsq,success)=leastsq(mglresiduals,coefs,args=(ydata,xdata,N),maxfev=50000)
    return [mglval(xdata,lsq,N),lsq]

def windowavg(data,N):
    return convolve(hanning(N)/sum(hanning(N)),data,mode='same')

def findpeakbounds(data,pkind,lowbnd,upbnd):
    fd=0.7
    [s,e]=[max(pkind-15,0),min(pkind+16,len(data))]
    mean=(max(data[s:e])-min(data[s:e]))*0.9+min(data[s:e])
    
#Left bound
    i=pkind
    #Step over plateau
    while fabs(data[i]-data[i-1]) < fd or data[i]>mean:
        if i>1: i-=1
        else:   break
    #Step down side
    while (data[i-1]-data[i]) < -1*fd:
        if i>1: i-=1
        else:   break
    lowbnd=max([pkind-50,lowbnd])
    start=max(i,lowbnd)

#Right bound
    i=pkind
    n=len(data)
    #Step over plateau
    while fabs(data[i]-data[i+1]) < fd or data[i]>mean:
        if i<n-2: i+=1
        else:     break
    #Step down side
    while (data[i+1]-data[i]) < -1*fd:   
        if i<n-2: i+=1
        else:     break
    upbnd=min([pkind+50,upbnd])
    end=min([i,upbnd])

    #check if not enough data was given using above algo, and artificially expand the bounds.
    d=end-start
    if d<15:
        end+=int((15-d)/2)+2
        start-=int((15-d)/2)-2

    return [start,end]

#find all 'clusters' of peaks
def findclusterbounds(data,pkinds):
    fdup=list()
    fddown=list()
    groups=list()
    fd=0.6
    i=25
    while i < len(data):
        found=False
        while (data[i]-data[i-1]) > fd:
            if found==False:
                found=True
                start=i
            i+=1
        if found==True:
            found=False
            for j in range(len(pkinds)):
                if abs(pkinds[j]-i)<10:
                    found=True
                    primepeak=j
                    break 
            if found==True:
                fdup.append(start)
                
                peaksowned=list()
                while i < len(data):
                    count=0
                    while i<len(data) and fabs(data[i]-data[i-1]) > fd:
                        count+=1
                        if primepeak < len(pkinds) and i>pkinds[primepeak]:
                            peaksowned.append(pkinds[primepeak])
                            primepeak+=1
                        i+=1
                    if count==0: #Plateau (step over)
                        i+=1
                        continue
                    count=0
                    while i<len(data) and fabs(data[i]-data[i-1]) < fd:
                        if count==0:
                            end=i
                        i+=1
                        count+=1
                    if count>10:
                        fddown.append(end)
                        if len(peaksowned)>1:
                            groups.append(peaksowned)
                        else:
                            fdup.pop()
                            fddown.pop()
                        break
                del peaksowned
        i+=1
        
    return (fdup,fddown,groups)

def extractpks(coefs,n):
    ilsq=list()
    for i in range(n):
        ilsq.append(list())
        for j in range(5):
            ilsq[i].append(coefs[i*5+j])
        ilsq[i].append(0)
        ilsq[i].append(0)
    return ilsq

#################
#given 3 points (xs,ys) to fit to, 
#returns coefficients a,b,c for a quadratic
def makequadratic(xs,ys): 
    mtx=[[x*x,x,1] for x in xs]
    return linalg.solve(mtx,ys)

#################
#given a ton of points (xs,ys) to fit to, 
#returns coefficients a,b,c for a quadratic
def makefitquadratic(xs,ys): 
    mtx=[[x*x,x,1] for x in xs]
    return linalg.lstsq(mtx,ys)[0]

def makefitlinear(xs,ys):
    mtx=[[x,1] for x in xs]
    return linalg.lstsq(mtx,ys)[0]

#################
#returns the y value of a quadratic function y=ax^2+bx+c
def quadval(x,(a,b,c)):
    return a*x*x+b*x+c

def lineval(x,(a,b)):
    return a*x+b

#################
#returns the y value of a gaussian function y=a*exp(-(b-x)^2/c)
def gaussianval(x,(a,b,c)):
    return a*exp(-(x-b)**2/c)

def gaussianres(coefs,y,x):
    return y-gaussianval(x,coefs)

def makefitgaussian(xs,ys,coefs):
    (lsq,success)=leastsq(gaussianres,coefs,args=(ys,xs),maxfev=10000)
    return lsq
