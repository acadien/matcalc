from itertools import imap
import cmath
import operator
#==================================================================
def dot(a,b):
    #a,b 3-vectors
    #returns the dotproduct of a and b.
    return sum(imap(operator.mul,a,b))

ci=-2.0*cmath.pi*complex(0,1) #Complex i
#==================================================================
def structurefactor(atoms,kvec):
    #atoms: list of atoms[N][3]
    #kvec: the k-vector to perform th
    return sum([cmath.exp(ci*dot(kvec,r)) for r in atoms])

#==================================================================
#Tools for Structure factor distribution calculations
ohtoone=[i/1000.0 for i in range(1001)] #includes 0 and 1
onetooh=[i/1000.0 for i in range(1000,-1,-1)] #includes 1 and 0
def SFdistro(atoms,a):
    #Returns the structure factor (intensity) for a standard set of
    #kvectors going from [0,0,0] to [0,0,a] to [0,a,a] to [0,0,0] to [a,a,a]
    
    SF=[structurefactor(atoms,[0*a,0*a,i*a]) for i in ohtoone]
    SF+=[structurefactor(atoms,[0*a,i*a,1.0*a]) for i in ohtoone]
    SF+=[structurefactor(atoms,[0*a,i*a,i*a]) for i in onetooh]
    SF+=[structurefactor(atoms,[i*a,i*a,i*a]) for i in ohtoone]
    return SF
