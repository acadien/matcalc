#!/usr/bin/python

#Writes an angular pair correlation to a file
def writeAngularPairCor(fname,header,minlen,maxlen,bins,vals):
    apcfile=open(fname,"w")
    apcfile.write(header.strip("\n")+"\n")
    apcfile.write(str(minlen)+" "+str(maxlen)+"\n")
    apcfile.write("\n".join([str(b)+" "+str(v) for b,v in zip(bins,vals)])+"\n")
    apcfile.close()
        
def readAngularPairCor(fname):
    apcdata=open(fname,"r").readlines()
    header=apcdata.pop(0)
    minlen,maxlen=map(float,apcdata.pop(0).split())
    bins,vals=zip(*[map(float,line.split()) for line in apcdata])
    return header,minlen,maxlen,bins,vals
