#!/usr/bin/python
import numpy as np
import glob,os

path="./"
for infile in glob.glob(os.path.join(path,'*.chi')):
    
    [Null,chifilename]=infile.split('/')
    
#chifilename="C-1_004.chi"

    outfilename=chifilename.split('.')[0]+".dat"

    f=open(chifilename,'r')
    i=0
    for line in f:
        if i<3:
            i+=1
            continue
        line=line.lstrip().rstrip()
        if i==3:
            i+=1
            length=int(line)
            break
    xxs=np.zeros(length)
    yys=np.zeros(length)
    i=0;
    for line in f:
        line=line.lstrip().rstrip()
        [xxs[i],NULL,yys[i]]=line.split(' ')
        i+=1
    f.close()

    f=open(outfilename,'w')
    data=zip(xxs,yys)
    [f.write(str(i[0])+" "+str(i[1])+"\n") for i in data]
    f.close()
