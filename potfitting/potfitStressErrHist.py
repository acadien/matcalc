#!/usr/bin/python

#mine
import plotRemote as pr

#notmine
import sys
import pylab as pl
from math import *

#Histograms the energy error results from potfit
#if the force database is given, different configurations can be labeled
#alter line labeled below to change how this is done

stressDirn=["xx","yy","zz","xy","yz","zx"]

def usage():
    print "%s <pf stress error file> <optional:force_DB> <optional: cutoff (=0.1)>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

sfile=open(sys.argv[1],"r").readlines()[2:]
fdb=-1
if len(sys.argv)>=3:
    fdb=["/".join(line.split("ascwork")[-1].split()[0].split("/")) for line in open(sys.argv[2],"r").readlines() if "ifconf" in line]

mxerr=0.05
if len(sys.argv)==4:
    mxerr=float(sys.argv[3])
sfile = map(lambda x: x.split(),sfile)
absdelS=map(lambda x: 0.0 if "inf" in x[5] else fabs(float(x[5])*float(x[4])),sfile)

if fdb!=-1:
    #Alter these lines to affect how labeling of configs is done
    eosDelS=[i for i,j in zip(absdelS,fdb) if "eos" in j]
    defectDelS=[i for i,j in zip(absdelS,fdb) if "defect" in j]
    meltDelS=[i for i,j in zip(absdelS,fdb) if "heat" in j or "cool" in j]
    feedDelS=[i for i,j in zip(absdelS,fdb) if "feedback" in j]
    strainDelS=[i for i,j in zip(absdelS,fdb) if "strain" in j]
    slowDelS=[i for i,j in zip(absdelS,fdb) if "slow" in j]
    print "EOSs: ",len(eosDelS)
    print "Melt/Quench: ",len(meltDelS)
    print "Feedback: ",len(feedDelS)
    print "Defects: ",len(defectDelS)
    print "Elastic: ",len(strainDelS)
    print "SlowQ: ",len(slowDelS)
    toHist = [i for i in [eosDelS,meltDelS,feedDelS,defectDelS,strainDelS,slowDelS] if len(i)>0]
    toHistL = [j for i,j in zip([eosDelS,meltDelS,feedDelS,defectDelS,strainDelS,slowDelS],["EOS","Melt","Feedback","Defects","Elastic","SlowQ"]) if len(i)>0]
    
    pl.hist(toHist,20,label=toHistL,histtype='barstacked',fill=True)
else:
    pl.hist(absdelS,20)
pl.xlabel("$|\Delta S|$")
pl.ylabel("count")
pl.legend(loc=0)

print "Worst fits: \ncnt   %-38.38s dirN     |delS|"%"config#"
cnt=0
for i,val in enumerate(reversed(sorted(absdelS))):
    if val<mxerr:
        break
    if fdb!=-1:
        ix=absdelS.index(val)
        print "%4d, %-40.40s %-3.3s | %f"%(ix,fdb[int(ix/6)],stressDirn[ix%6],absdelS[ix])
    else:
        print i," ",absdelS.index(val)/6," \t",val
    cnt=i
print "\n%d configurations with error above %f"%(cnt,mxerr)

pr.prshow("stressErrorHist.png")
#pl.savefig("/home/acadien/Dropbox/GeData/hist_energy_error.png")

