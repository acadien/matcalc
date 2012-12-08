#!/usr/bin/python

def write(fname,header,bins,vals):
    rdffile=open(fname,"w")
    rdffile.write(header.strip("\n")+"\n")
    rdffile.write("\n".join([str(b)+" "+str(v) for b,v in zip(bins,vals)])+"\n")
    rdffile.close()
        
def read(fname):
    rdfdata=open(fname,"r").readlines()
    header=rdfdata.pop(0)
    bins,vals=zip(*[map(float,line.split()) for line in rdfdata])
    return header,bins,vals
