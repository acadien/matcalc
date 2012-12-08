#!/usr/bin/python

#Writes the coordination information to a file
def write(fname,header,mncn,mxcn,labels,avgs,cnhists):
    seperator="|||"
    cnfile=open(fname,"w")
    cnfile.write(header.strip("\n")+"\n")
    cnfile.write(str(mncn)+" "+str(mxcn)+"\n")
    cnfile.write(seperator.join(labels)+"\n")
    cnfile.write(" ".join(map(str,avgs))+"\n")
    for i in cnhists:
        cnfile.write(" ".join(map(str,i))+"\n")
    cnfile.close()

def read(fname):
    seperator="|||"
    cndata=open(fname,"r").readlines()
    header=cndata.pop(0)
    mncn,mxcn=map(int,cndata.pop(0).split())
    labels=cndata.pop(0).split(seperator)
    avgs=map(float,cndata.pop(0).split())
    cnhists=map(lambda x: map(float,x.split()),cndata)
    return header,mncn,mxcn,labels,avgs,cnhists
    
