#!/usr/bin/python

def writePairCor(fname,header,bins,vals):
    apcfile=open(fname,"w")
    apcfile.write(header.strip("\n")+"\n")
    apcfile.write("\n".join([str(b)+" "+str(v) for b,v in zip(bins,vals)])+"\n")
    apcfile.close()
        
def readPairCor(fname):
    apcdata=open(fname,"r").readlines()
    header=apcdata.pop(0)
    bins,vals=zip(*[map(float,line.split()) for line in apcdata])
    return header,bins,vals
