#!/usr/bin/python

from numpy import matrix,linalg,array

#Given the poscar (opened and read)
#Returns: basis,atypes,atoms,head,poscar
def read(poscar,frac_coord=False,selDynFlag=False):
    if type(poscar)!=type(list()):
        poscar=open(poscar,"r").readlines()

    if len(poscar)<3:
        print "Error\nSomething wrong with poscar or not passed properly to read poscar."
        return [-1]*5

    head=poscar.pop(0)
    scale=float(poscar.pop(0))

    v1=[float(i)*scale for i in poscar.pop(0).split()]
    v2=[float(i)*scale for i in poscar.pop(0).split()]
    v3=[float(i)*scale for i in poscar.pop(0).split()]
    try:
        vals=poscar.pop(0).split()
        atypes=[int(i) for i in vals]
    except ValueError:
        elements=vals
        atypes=[int(i) for i in poscar.pop(0).split()]
    N=sum(atypes)

    line = poscar.pop(0)
    if "Sel" in line or "sel" in line: #selective dynamics
        poscar.pop(0)

    seldyn=list()
    if selDynFlag:
        ax,ay,az,sd1,sd2,sd3=map(list,zip(*[line.split() for line in poscar[:N]]))        
        ax=map(float,ax)
        ay=map(float,ay)
        az=map(float,az)
        seldyn=zip(sd1,sd2,sd3)
    else:
        ax,ay,az=map(list,zip(*[map(float,line.split()[:3]) for line in poscar[:N]]))        

    if not frac_coord:
        for i in range(len(ax)):
            a,b,c = array(ax[i]), array(ay[i]), array(az[i])
            ax[i]=v1[0]*a+v2[0]*b+v3[0]*c
            ay[i]=v1[1]*a+v2[1]*b+v3[1]*c
            az[i]=v1[2]*a+v2[2]*b+v3[2]*c

        center=(sum(ax)/len(ax)+sum(ay)/len(ay)+sum(az)/len(az))/3
        if center<1: #in fractional coordinates, convert to cartesian
            print "WARNING: poscarIO.read() did not convert to fractional coordinates, something probably wrong with POSCAR."
    
    atoms=zip(ax,ay,az)
    basis=[v1,v2,v3]
    poscar=poscar[N:]    
    
    if selDynFlag:
            return basis,atypes,atoms,head,poscar,seldyn
    return basis,atypes,atoms,head,poscar

#Takes 1 valid poscar from the input and returns it with the input poscar
def split(poscarin):
    if len(poscarin)<3:
        return [-1]*2

    poscarout=poscarin[:7]
    Natoms=sum([int(i) for i in poscarin[5].split()])
    poscarin=poscarin[7:]
    for i in range(Natoms):
        poscarout.append(poscarin.pop(0))

    return poscarin,poscarout


def write(poscarName,basis,atoms,types,header,seldyn=None,ratio=1.0,frac=True,autoFrac=True,velocities=None):

    if len(types)==len(atoms):
        #fix types arrangement, reorder atoms as necessary
        atoms,types = zip(*sorted(zip(atoms,types),key=lambda x:x[1]))
        sm=min(types)
        bg=max(types)
        ntypes=[types.count(i) for i in range(sm,bg+1)]
    else:
        ntypes=types

    if seldyn==None:
        seldyn=[""]*len(atoms)
        seldynTag=""
    else:
        seldyn=[" ".join(i) for i in seldyn]
        seldynTag="Selective Dynamics\n"

    #If atoms are not in fractional coordinates, convert them
    if autoFrac==True:
        com = array(atoms).sum().sum()/3/len(atoms)
        if com > 1.0 or com < 0.0 :
            frac=False
        else:
            frac=True

    #Convert if necessary (re-using the frac flag)
    if not(frac):
        A = matrix(basis)
        atoms=[linalg.solve(A.T,array(atom))[:] for atom in atoms]

    #Make the POSCAR
    data=""
    data+=" ".join(header.split("\n"))+"\n"
    data+=" % 16.16f\n"%ratio
    for vec in basis:
        data+=" % 5.12e % 5.12e % 5.12e\n"%(vec[0],vec[1],vec[2])
    data+=" "+" ".join([str(i) for i in ntypes])+"\n"
    data+=seldynTag
    data+="Direct\n"
    for sd,atom in zip(seldyn,atoms):
        data+=" % 5.12e % 5.12e % 5.12e %s\n"%(atom[0],atom[1],atom[2],sd)

    #Velcity
    if velocities!=None:
        print "here",velocities,velocities==None
        data+="\n"
        for v in velocities:
            data+=" % 5.12e % 5.12e % 5.12e\n"%(v[0],v[1],v[2])

    poscar=open(poscarName,"w")
    poscar.write(data)
    return poscar
