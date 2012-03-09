#!/usr/bin/python

import os

def dir(remotedir,localdir):
    #os.system("scp -r 
    L = ['scp', '-r','gmice.gmu.edu:/home/acadien/%s/'%remotedir, localdir]
    os.spawnvpe(os.P_WAIT, 'scp', L, os.environ)


def file(remotedir,localdir,localfil):
    #os.system("scp gmice.gmu.edu:/home/acadien/%s/%s %s"%(remotedir,localfil,localdir))  
    L = ['scp','gmice.gmu.edu:/home/acadien/%s/%s'%(remotedir,localfil), localdir]
    os.spawnvpe(os.P_WAIT, 'scp', L, os.environ)

