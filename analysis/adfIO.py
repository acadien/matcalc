#!/usr/bin/python

#Writes an angular distribution function to a file
def write(fname,header,minlen,maxlen,bins,vals):
    adffile=open(fname,"w")
    adffile.write(header.strip("\n")+"\n")
    adffile.write(str(minlen)+" "+str(maxlen)+"\n")
    adffile.write("\n".join([str(b)+" "+str(v) for b,v in zip(bins,vals)])+"\n")
    adffile.close()
        
def read(fname):
    adfdata=open(fname,"r").readlines()
    header=adfdata.pop(0)
    minlen,maxlen=map(float,adfdata.pop(0).split())
    bins,vals=zip(*[map(float,line.split()) for line in adfdata])
    return header,minlen,maxlen,bins,vals
