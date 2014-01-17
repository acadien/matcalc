#!/usr/bin/python

#mine
import plotRemote as pr

#notmine
import sys
import pylab as pl

#Histograms the energy error results from potfit
#if the force database is given, different configurations can be labeled
#alter line labeled below to change how this is done

def usage():
    print "%s <pf energy error file> <optional:force_DB> <optional: cutoff (=0.1)>"%sys.argv[0].split("/")[-1]

if len(sys.argv)<2:
    usage()
    exit(0)

efile=open(sys.argv[1],"r").readlines()[2:]
fdb=-1
if len(sys.argv)>=3:
    fdb=["/".join(line.split("home")[1].split("OUTCAR")[0].split("/")[6:]) for line in open(sys.argv[2],"r").readlines() if "ifconf" in line]

mxerr=0.1
if len(sys.argv)==4:
    mxerr=float(sys.argv[3])

absdelE=map(lambda x:float(x.split()[5]),efile)
if fdb!=-1:
    #Alter these lines to affect how labeling of configs is done
    eosDelE=[i for i,j in zip(absdelE,fdb) if "eos" in j]
    defectDelE=[i for i,j in zip(absdelE,fdb) if "defect" in j]
    meltDelE=[i for i,j in zip(absdelE,fdb) if "heat" in j or "cool" in j]
    feedDelE=[i for i,j in zip(absdelE,fdb) if "feedback" in j]
    strainDelE=[i for i,j in zip(absdelE,fdb) if "strain" in j]
    print "EOSs: ",len(eosDelE)
    print "Defects: ",len(defectDelE)
    print "Melt/Quench: ",len(meltDelE)
    print "Feedback: ",len(feedDelE)
    print "Elastic: ",len(strainDelE)

    labels=["EOS","Melt","Feedback","Defects","Elastic"]
    toHist = [i for i in [eosDelE,meltDelE,feedDelE,defectDelE,strainDelE] if len(i)>0]
    toLabel = [labels[j] for j,i in enumerate([eosDelE,meltDelE,feedDelE,defectDelE,strainDelE]) if len(i)>0]
    pl.hist(toHist,20,label=toLabel,histtype='barstacked',fill=True)
else:
    pl.hist(absdelE,20)
pl.xlabel("$|\Delta E|$")
pl.ylabel("count")
pl.legend(loc=0)

print "Worst fits: \ncnt   %-48.48s |delE|"%"config#"
cnt=0
for i,val in enumerate(reversed(sorted(absdelE))):
    if val<mxerr:
        break
    if fdb!=-1:
        print "%4d, %-50.50s, %f"%(absdelE.index(val),fdb[absdelE.index(val)],val)
    else:
        print i," ",absdelE.index(val)," \t",val
    cnt=i
print "\n%d configurations with error above %f"%(cnt,mxerr)

pr.prshow("energyErrorHist.png")
#pl.savefig("/home/acadien/Dropbox/GeData/hist_energy_error.png")

