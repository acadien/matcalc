#!/usr/bin/python
#Calculation tool for POSCARS and LAMMPS dumps.

import sys
#mine
import plotRemote as pr

from colors import float2rgb
from scipy import array
from numpy import linspace,savetxt
import matplotlib.cm as cm
import pylab as pl
from math import sqrt
import argparse
import re
#mine
from orderParam import *
import poscarIO
import outcarIO
import lammpsIO

orderParams={"CN":  coordinationNumber, 
             "BO":  bondOrientation , 
             "RDF": radialDistribution , 
             "ADF": angleDistribution , 
             "BA":  bondAngleCorr , 
             "SF":  structureFactor, 
             "TET": tetrahedral, 
             "TN":  translational}

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''\
Apply some order parameter calculation to a POSCAR, OUTCAR or LAMMPS Dump, then attempt to plot it or print out the result.\n
==========================================================
= Order Parameter must be one of:                        =
=   CN  : Coordination Number                            =
=   BO  : Bond Orientation (Q_l), requires -lval be set  =
=   RDF : Radial Distribution Function                   =
=   ADF : Angular Distribution Function                  =
=   BA  : Bond Angle Correlation, requires -gval be set  =
=   SF  : Structure Factor                               =
=   TET : Tetrahedral                                    =
=   TN  : Translational (tao)                            = 
==========================================================
''',epilog='''\
Examples:
charter.py RDF POSCAR1 POSCAR2 POSCAR2 -A
charter.py ADF OUTCAR -N 1-200 -rcut 3.5
charter.py CN dumpfile.dat -N 25,36 -rcut 6.5
''')

#Need to implement test of file to determine if it is a lammps dump or poscar chart
#Also implement for OUTCAR files...

parser.add_argument('op',choices=orderParams.keys(),help="Order Parameter selection.")
parser.add_argument('fileNames',nargs='*',help="A list of POSCAR files (space delim)")
parser.add_argument('-rcut',help='Radius Cutoff distance when selecting neighbors, if not provided the first shell minimum is used.',type=float)
parser.add_argument('-lval',dest='lval',help='l-value, used in BO (Ql) calculation',type=int)
parser.add_argument('-gval',dest='lval',help='g-value, used in BA (g) calculation',type=int)
parser.add_argument('-nosave','-S',dest='saveFlag',action='store_false',help='turn off saving to file',default='true')
parser.add_argument('-noplot','-P',dest='plotFlag',action='store_false',help='turn off plotting',default='true')
parser.add_argument('-avg','-A',dest='averageFlag',action='store_true',help='average POSCAR results if possible')
parser.add_argument('-debug','-D',dest='debug',action='store_true',help='turn on debugging of neighbor selection')
parser.add_argument('-N',dest='cfgNums',help='Configuration number list from dump/outcar, can be comma separated list, or dash separated, #1-#2, to get range (default -1)',type=str,default="-1")
#Need to implement this
parser.add_argument('-smooth','-sm',dest='windowSize',help='Turns on windowed averaging (smoothing) over all data sets',type=int,default=10)

args= parser.parse_args()
op = args.op
fileNames = args.fileNames

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
        
lval = args.lval
cfgNums = args.cfgNums
if "-" in cfgNums and len(cfgNums)>2:
    a,b=map(int,cfgNums.split("-"))
    cfgNums=range(a,b)
else:
    cfgNums = map(int,args.cfgNums.split(","))


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

#======================================================
#       Gather the order parameters
#======================================================
neighbs=None #maybe do something with neighbors option later
orderVals=list()
for pn in fileNames:

    if "POSCAR" in pn or "CONTCAR" in fileNames[0]:
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
          "ADF":[r"$\theta(deg)$",               "#_Bond_Angles"],
          "BA" :[r"R$(\AA)$",                    r"$G_{%s}(r)$"%str(lval)],
          "SF" :[r"Q$(\AA^{-1})$",               r"S(Q)"]}

if args.averageFlag:
    if op not in ["BO","CN","TN","TET"]:
        t = zip(*[orderVals[i][1] for i in range(len(orderVals))])
        avgy = map(lambda x:sum(x)/len(x),t)
        avgx = orderVals[0][0]
    elif op=="CN":
        hd,rcuts=zip(*orderVals)
        hd =  [i for i in flatten(list(hd))] #flatten for averaging
        avgy,avgx,dummy=pl.hist(hd,bins=range(0,16),visible=False,normed=True)
        avgx = avgx[:-1]

    pl.xlabel(xylabels["CN"][0])
    pl.ylabel(xylabels["CN"][1])
    if args.saveFlag:
        prefix=""
        if sameDir:
            prefix=fileDirs[0]

        if op in ["TN","TET"]:
            savetxt(prefix+"AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ")
        elif op in ["ADF"]:
            savetxt(prefix+"AVERAGE."+op+str(args.rcut),array([avgx,avgy]).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
        else:
            savetxt(prefix+"AVERAGE."+op+str(lval),array([avgx,avgy]).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
    if args.plotFlag:
        pl.plot(avgx,avgy)

elif args.saveFlag:
    for ov,pn in zip(orderVals,fileNames):
        if op in ["TN","TET"]:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ")

        elif op in ["ADF"]:
            if args.rcut==None:
                savetxt(pn+"."+op,array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")
            else:
                savetxt(pn+"."+op+str(args.rcut),array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")

        elif op in ["BO"]:
            vals,bins,dummy = pl.hist(ov,bins=int(sqrt(len(ov)))*5,normed=True,visible=False)
            savetxt(pn+"."+op+str(lval),array(zip(bins[:-1],vals)),delimiter=" ",header=" ".join(xylabels[op]),comments="")

        else:
            savetxt(pn+"."+op+str(lval),array(ov).T,delimiter=" ",header=" ".join(xylabels[op]),comments="")

if not args.plotFlag:
    exit(0)

if op not in ["BO","CN","TN","TET"] and not args.averageFlag:
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
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)/len(ov)

elif op=="TN":
    print "Average Translational Order <tao> ="
    for i,ov in enumerate(orderVals):
        print fileNames[i],"\t\t",sum(ov)

if args.averageFlag:
    pl.legend(["Average "+" ".join(fileNames)],loc=0,fontsize=10)
else:
    if len(labels)>0:
        pl.legend(labels,loc=0,fontsize=10)
    else:
        pl.legend(fileNames,loc=0,fontsize=10)

#Don't plot for some values
if op not in ["TN","TET"]:
    pr.prshow("%s_chart_%s.png"%(sys.argv[1],op))
