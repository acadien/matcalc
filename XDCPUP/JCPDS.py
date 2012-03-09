#!/usr/bin/python

#A class for holding JCPDS information

class JCPDS:
    
    def __init__(self):
        self.comment=""
        
        #Thermodynamics
        self.bulk=0.0
        self.bulkP=0.0

        #Unit cell parameters
        self.Symmetry=""
        self.A=0.0
        self.B=0.0
        self.C=0.0
        self.alpha=0.0
        self.beta=0.0
        self.gamma=0.0
        
        #Dspacings
        self.hkls=list()
        self.dspacings=list()
        self.intensities=list()
        print "init_jcpds"

#Given the JCPDS File, reads in the d-spacings, intensities and cell parameters
#Returns a JCPDS object with all of the corresponding information.
    def load_jcpds(self,fname):
        print "load_jcpds"
        jcpds_lines=open(fname).readlines()
    
        for line in jcpds_lines:
            words=line.split()
        
            if words[0]=="COMMENT:":
                self.comment=line
            elif words[0]=="K0:":
                self.bulk=float(words[1])
            elif words[0]=="K0P:":
                self.bulkp=
            elif words[0]=="SYMMETRY:":
