#!/usr/bin/python
#Calculation tool for POSCARS and LAMMPS dumps.

import sys
#mine - must be imported before scipy/numpy/mpl
import plotRemote as pr

from colors import float2rgb,vizSpec
from scipy import array
from numpy import linspace,savetxt,histogram
import matplotlib.cm as cm
import pylab as pl
from math import sqrt
import argparse
import re
import matplotlib
#3D
from mpl_toolkits.mplot3d import axes3d
#mine
from orderParam import *
import poscarIO
import outcarIO
import lammpsIO

orderParams={"CN":  coordinationNumber, 
             "BO":  bondOrientation , 
             "RDF": radialDistribution , 
             "FRDF": radialDistribution ,
             "ADF": angleDistribution , 
             "RAF": radangDistribution ,
             "BA":  bondAngleCorr , 
             "SF":  structureFactor, 
             "SF0": structureFactor0,
             "TET": tetrahedral, 
             "TN":  translational}

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''\
Apply some order parameter calculation to a POSCAR, OUTCAR or LAMMPS Dump, then attempt to plot it or print out the result.\n
==========================================================
= Order Parameter must be one of:                        =
=   CN  : Coordination Number                            =
=   BO  : Bond Orientation (Q_l), requires -lval be set  =
=   RDF : Radial Distribution Function                   =
=   FRDF: Fuzzy Radial Distribution Function             =
=   ADF : Angular Distribution Function                  =
=   RAF : Radial distribution by bond angle              =
=   BA  : Bond Angle Correlation, requires -gval be set  =
=   SF  : Structure Factor                               =
=   SF0 : Structure Factor (k=0)                         =
=   TET : Tetrahedral                                    =
=   TN  : Translational (tao)                            = 
==========================================================
''',epilog='''\
Examples:
charter.py RDF POSCAR1 POSCAR2 POSCAR2 -A
charter.py ADF OUTCAR -N 1-200 -rcut 3.5
charter.py CN dumpfile.dat -N 25,36 -rcut 6.5
''')

parser.add_argument('op',choices=orderParams.keys(),help="Order Parameter selection.")
parser.add_argument('fileNames',nargs='*',help="A list of POSCAR files (space delim)")
parser.add_argument('-rcut',help='Radius Cutoff distance when selecting neighbors, if not provided the first shell minimum is used.',type=float)
parser.add_argument('-lval',dest='lval',help='l-value, used in BO (Ql) calculation',type=int)
parser.add_argument('-gval',dest='lval',help='g-value, used in BA (g) calculation',type=int)
parser.add_argument('-nosave','-S',dest='saveFlag',action='store_false',help='turn off saving to file',default='true')
parser.add_argument('-noplot','-P',dest='plotFlag',action='store_false',help='turn off plotting',default='true')
parser.add_argument('-avg','-A',dest='averageFlag',action='store_true',help='average results if possible')
parser.add_argument('-debug','-D',dest='debug',action='store_true',help='turn on debugging of neighbor selection')
parser.add_argument('-N',dest='cfgNums',help='Configuration number list from dump/outcar, can be comma separated list, or dash separated, #1-#2, to get range (default -1), can also be \'all\'',type=str,default="-1")
parser.add_argument('-saveas',dest='saveFileName',help='Filename under which to save the data to',type=str,default='None')
parser.add_argument('-nbins',dest="nbins",help="Number of bins to use in the operation",type=int,default=None)
parser.add_argument('-evolve',dest='nEvolve',help='Track the time evolution of a parameter, average over nEvolve stages',type=int,default=None)
parser.add_argument('-stagger',dest='stagger',action='store_true',help='staggers results when plotting, useful with -evolve setting',default=False)
parser.add_argument('-neighbors',dest='neighbFile',help='A file containing a nieghbor list',type=str,default=None)
args= parser.parse_args()

op = args.op
fileNames = args.fileNames
saveFileName = args.saveFileName
nEvolve = args.nEvolve
stagger = args.stagger
neighbFile = args.neighbFile
lval = args.lval
cfgNums = args.cfgNums
nbins = args.nbins

if op=="RDF":
    lval=nbins

#grab the directory for each file
fileDirs = list()
for fn in fileNames:
    fns = fn.split("/")
    if len(fns)==1:
        fileDirs.append("./")
    else:
        fileDirs.append("/".join(fns[:-1])+"/")
#Are all the file in the same directory?
sameDir = False
if fileDirs.count(fileDirs[0]) == len(fileDirs):
    sameDir = True

#Loops in cfgNums
if "-" in cfgNums and len(cfgNums)>2:
    a,b=map(int,cfgNums.split("-"))
    cfgNums=range(a,b)
elif "all" in cfgNums or "All" in cfgNums:
    pass
else:
    cfgNums = map(int,args.cfgNums.split(","))

if nEvolve!=None and args.cfgNums=="-1":
    print "Must specify configurations using -N in order to use Evolve command"
    exit(0)

if op == "BA" and args.lval == None:
    print "Error: g-value (-gval) must be set for BA order parameter"
    exit(0)

if op == "BO" and args.lval == None:
    print "Error: l-value (-lval) must be set for BO order parameter"
    exit(0)
    
#Sort POSCAR names only based on the numbers in them
try:
    fileNumbers=map(lambda x:float("".join(re.findall('\d+',x))),fileNames)
    if len(fileNumbers) == len(fileNames):
        fileNames=zip(*sorted(zip(fileNames,fileNumbers),key=lambda x:x[1]))[0]
except ValueError:
    pass

#Handle Neighbor Files
neighbs=None 
if neighbFile != None:
    neighbs=neighborIO.read(neighbFile)

#======================================================
#       Gather the order parameters
#======================================================
orderVals=list()
for pn in fileNames:

    if "POSCAR" in pn or "CONTCAR" in pn:#fileNames[0]:
        poscar=open(pn,"r").readlines()
        [basis,atypes,atoms,head,poscar] = poscarIO.read(poscar)
        orderVals.append(orderParams[op]( \
                array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))

    elif "OUTCAR" in pn:
        if len(cfgNums)==1:
            TE,stress,basis,atoms,forces,types = outcarIO.outcarReadConfig(pn,cfgNums)
            orderVals.append( orderParams[op]( \
                    array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))
        else:
            TEs,stresss,basiss,atomss,forcess,typess = outcarIO.outcarReadConfig(pn,cfgNums)
            for basis,atoms in zip(basiss,atomss):
                orderVals.append(orderParams[op]( \
                    array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))

    else: #LAMMPS DUMP
        if len(cfgNums)==1:
            basis,types,atoms =  lammpsIO.readConfig(pn,cfgNums)
            orderVals.append( orderParams[op]( \
                    array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))
        else:
            basiss,types,atomss = lammpsIO.readConfig(pn,cfgNums)
            for basis,atoms in zip(basiss,atomss):
                orderVals.append(orderParams[op]( \
                    array(atoms),array(basis),l=lval,neighbs=neighbs,rcut=args.rcut,debug=args.debug))

#======================================================
#                       Plot!
#======================================================
labels=[]
if lval==None:
    lval=""
xylabels={"BO" :[r"BondOrder($Q_%s$)"%str(lval), r"$P(Q_%s)$"%str(lval)],
          "CN" :[r"Coordination_Number",         r"P(CN)"],
          "RDF":[r"R$(\AA)$",                    "#_Bonds"],
          "FRDF":[r"R$(\AA)$",                    "#_Bonds"],
          "ADF":[r"$\theta(deg)$",               "#_Bond_Angles"],
          "BA" :[r"R$(\AA)$",                    r"$G_{%s}(r)$"%str(lval)],
          "SF" :[r"Q$(\AA^{-1})$",               r"S(Q)"],
          "RAF":[r"temp",r"temp2"],
          "TET":[r"$q_t$","count"]}

#Special case for handling 3D RAF data.
if op == "RAF":
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')    
    if args.averageFlag:
        N=float(len(orderVals))
        for k,((adfVals,rdfVals),bins) in enumerate(orderVals):
            X=[[i]*len(rdfVals) for i in adfVals]
            Y=[rdfVals for i in range(len(adfVals))]
            if k==0:
                Z=array([[b/N for b in a]for a in bins])
            else:
                Z+=array([[b/N for b in a]for a in bins])
        Z/=len(orderVals)
    else:
        for k,((adfVals,rdfVals),bins) in enumerate(orderVals):
            X=[[i]*len(rdfVals) for i in adfVals]
            Y=[rdfVals for i in range(len(adfVals))]
            Z=bins

    ax.plot_wireframe(X,Y,Z,color='black',lw=0.5)
    ax.contour(X, Y, Z,10,cmap=cm.coolwarm,linewidths=[5]*10)
    pl.xlabel("theta")
    pl.ylabel("r")
    pl.show()
    exit(0)

#Special Case for FRDF histogrammic RDF
if op == "FRDF":
    from numpy import bincount,array,zeros
    import numpy as np
    from datatools import wsmooth

    t = zip(*[[j*100 for j in orderVals[i][1]] for i in range(len(orderVals))])
    avgx = orderVals[0][0]
    a=zeros([len(t),2000])
    mx=0
    for i,v in enumerate(t):
        a[i]=bincount(v,minlength=2000)
        a[i]=wsmooth(a[i],window_len=75)[0:2000]
        m=max(a[i])
        a[i]=np.log([j/m for j in a[i]])

    pl.imshow(a.T,origin='lower',cmap=pl.get_cmap("Blues"))
    pl.yticks(range(0,1000,100),range(10))
    pl.xticks(range(0,1001,100),range(11))

    pl.show()
    exit(0)

def savetxtWrapper(defaultFileName,data,delimiter=" ",header=" ",comments=""):
    try:
        if saveFileName is "None":
            savetxt(defaultFileName,data,delimiter=delimiter,header=header,comments=comments)
        else:
            savetxt(saveFileName,data,delimiter=delimiter,header=header,comments=comments)
    except TypeError:
        print "unable to save... moving along"


if nEvolve != None:
    [xs,ys]=zip(*orderVals)
    xs=xs[0]
    nSamples=len(ys)

    nDel=nSamples/nEvolve
    yEvs=list()
    lablel=list()
    for i in range(nEvolve):
        start, finish = i*nDel, (i+1)*nDel
        labels.append(str(start)+" - "+str(finish))
        yEvs.append(map(lambda x: sum(x)/len(x),zip(*ys[start:finish])))
        
    if not(stagger):
        for i in range(nEvolve):
            pl.plot(xs,yEvs[i],color=vizSpec(i/float(nEvolve)))
    else:
        yDel=( max(map(max,yEvs)) - min(map(min,yEvs)) )/2.
        for i in range(nEvolve):
            yvals=[y+yDel*i for y in yEvs[i]]
            pl.plot(xs,yvals,color=vizSpec(i/float(nEvolve)))
            pl.text((max(xs)-min(xs))*0.8,yDel*(i+0.2),labels[i])
        pl.yticks([])
    pl.xlabel(xylabels[op][0])
    pl.ylabel(xylabels[op][1])
    pr.prshow("%s_chart_evolve.png"%op)
    exit(0)

if args.averageFlag:
    if op not in ["BO","CN","TN","TET"]:
        t = zip(*[orderVals[i][1] for i in range(len(orderVals))])
        avgy = map(lambda x:sum(x)/len(x),t)
        avgx = orderVals[0][0]
    elif op in ["CN","TET"]:
        hd,rcuts=zip(*orderVals)
        hd =  [i for i in flatten(list(hd))] #flatten for averaging
        mn = min(hd)
        mx = max(hd)
        d = mx-mn
        mn -= 0.2*d
        mx += 0.2*d
        avgy,avgx,dummy=pl.hist(hd,bins=16,range=(mn,mx),visible=False,normed=True)
        avgx = avgx[:-1]

    pl.xlabel(xylabels[op][0])
    pl.ylabel(xylabels[op][1])
    if args.saveFlag:
        prefix=""
        if sameDir:
            prefix=fileDirs[0]

        if op in ["TN"]:
            savetxtWrapper(prefix+"AVERAGE."+op+str(lval),array([avgx,avgy]).T)
        elif op in ["ADF"]:
            savetxtWrapper(prefix+"AVERAGE."+op+str(args.rcut),array([avgx,avgy]).T,header=" ".join(xylabels[op]))
        else:
            savetxtWrapper(prefix+"AVERAGE."+op+str(lval),array([avgx,avgy]).T,header=" ".join(xylabels[op]))
    if args.plotFlag:
        pl.plot(avgx,avgy)

elif args.saveFlag:
    for ov,pn in zip(orderVals,fileNames):
        if op in ["TN","TET"]:
            vals,bins,dummy = pl.hist(ov[0],bins=int(sqrt(len(ov[0])))*2,normed=True,visible=False)
            savetxtWrapper(pn+"."+op+str(lval),array(zip(bins[:-1],vals)))

        elif op in ["ADF"]:
            if args.rcut==None:
                savetxtWrapper(pn+"."+op,array(ov).T,header=" ".join(xylabels[op]))
            else:
                savetxtWrapper(pn+"."+op+str(args.rcut),array(ov).T,header=" ".join(xylabels[op]))

        elif op in ["BO"]:
            vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,normed=True,visible=False)
            savetxtWrapper(pn+"."+op+str(lval),array(zip(bins[:-1],vals)),header=" ".join(xylabels[op]))
        elif op == "SF0":
            print "still figuring out how to do this calculation :("
        else:
            savetxtWrapper(pn+"."+op+str(lval),array(ov).T,header=" ".join(xylabels[op]))

if not args.plotFlag:
    exit(0)

if op not in ["BO","CN","TN","TET","SF0"] and not args.averageFlag:
    for ov in orderVals:
        pl.plot(ov[0],ov[1])
    pl.xlabel(xylabels[op][0])
    pl.ylabel(xylabels[op][1])

if op=="BO":
    print "Average Bond Order <BO%d> ="%lval
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)/len(ov)
        vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,visible=False)#,histtype='step')
        pl.plot(array(bins[:-1]),array(vals),label=fileNames[i])
    pl.xlabel(xylabels["BO"][0])
    pl.ylabel(xylabels["BO"][1])

elif op=="CN":
    hd,rcuts=zip(*orderVals)
    pl.hist(hd,bins=range(0,16),normed=True,histtype='bar',align='left',rwidth=0.8)
    pl.xticks(range(min(map(min,hd)),max(map(max,hd))+1))
    #for i,ov in enumerate(orderVals):
    #    cn=float(sum(ov[0]))/len(ov[0])
    #    labels.append("CN%3.3f rcut%3.3f %s"%(cn,rcuts[i],fileNames[i]))
    pl.xlabel(xylabels["CN"][0])
    pl.ylabel(xylabels["CN"][1])

elif op=="TET":
    print "Average Tetrahedral Order <Sg> ="
    avgTets=list()
    for i,ov in enumerate(orderVals):
        tet=ov[0]
        if len(orderVals)>1:
            avgTets.append(sum(tet)/len(tet))
        else:
            avgTets.append(tet)
        print sum(tet)/len(tet)
    print "min, max, avg"
    if len(orderVals)>1:
        print min(avgTets),max(avgTets),sum(avgTets)/len(avgTets)
    else:
        print min(avgTets[0]),max(avgTets[0]),sum(avgTets[0])/len(avgTets[0])
    pl.hist([ov[0] for ov in orderVals])

elif op=="TN":
    print "Average Translational Order <tao> ="
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)

elif op=="SF0":
    print "shit still figuring this out"
    exit(0)

if args.averageFlag:
    pl.legend(["Average "+" ".join(fileNames)],loc=0,fontsize=10)
else:
    if len(labels)>0:
        pl.legend(labels,loc=0,fontsize=10)
    else:
        pl.legend(fileNames,loc=0,fontsize=10)

#Don't plot for some values
if op not in ["TN"]:
    pr.prshow("%s_chart_%s.png"%(sys.argv[1],op))
