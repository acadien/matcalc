#!/usr/bin/python

#Read and write neighbor files

def write(fname,neighbs):
    nfile=open(fname,"w")
    nfile.writelines(map(lambda x: " ".join(map(str,x))+"\n",neighbs))

def read(fname):
    neighbs=list()
    for line in open(fname,"r").readlines():
        if line[0]=="#" or (len(neighbs)==0 and len(line)==0):
            continue
        
        if len(line.split())==0:
            neighbs.append([0])
        else:
            neighbs.append(map(int,line.split()))
    return neighbs
