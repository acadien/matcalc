#!/usr/bin/python

from numpy import matrix,linalg

def read(poscar,frac_coord=False):
    if len(poscar)<3:
        print "Error\nSomething wrong with poscar or not passed properly to read poscar."
        return [-1]*9

    head=poscar.pop(0)
    scale=float(poscar.pop(0))**(1./3)

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
    poscar.pop(0)

    ax=list()
    ay=list()
    az=list()

    if frac_coord:
        for line in poscar[:N]:
            [x,y,z]=[float(i) for i in line.split()]
            ax.append(x)
            ay.append(y)
            az.append(z)
    else:
        for line in poscar[:N]:
            [x,y,z]=[float(i) for i in line.split()]
            x1=v1[0]*x+v2[0]*y+v3[0]*z
            y1=v1[1]*x+v2[1]*y+v3[1]*z
            z1=v1[2]*x+v2[2]*y+v3[2]*z
            ax.append(x1)
            ay.append(y1)
            az.append(z1)
    poscar=poscar[N:]    
        
    return v1,v2,v3,atypes,ax,ay,az,head,poscar

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


def write(poscarName,basis,atoms,types,header,frac=True,ratio=1.0):

    if len(types)==len(atoms):
        #fix types arrangement, reorder atoms as necessary
        atoms,types = zip(*sorted(zip(atoms,types),key=lambda x:x[1]))
        sm=min(types)
        bg=max(types)
        ntypes=[types.count(i) for i in range(sm,bg+1)]

    #If atoms are not in fractional coordinates, convert them
    if frac:
        #Check if they're already fractional
        for d,atom in enumerate(atoms):
            if len([1 for i in atom if i>1.0 or i<0.0])>0:
                frac=False
                break
        #Convert if necessary (re-using the frac flag)
        if not(frac):
            A = matrix(basis)
            atoms=[linalg.solve(A,atom)[:] for atom in atoms]

    #Make the POSCAR
    data=""
    data+=" ".join(header.split("\n"))+"\n"
    data+="%3.3f\n"%ratio
    for vec in basis:
        data+="% 5.12e % 5.12e % 5.12e\n"%(vec[0],vec[1],vec[2])
    data+=" ".join([str(i) for i in types])+"\n"
    data+="Direct\n"
    for atom in atoms:
        data+="% 5.12e % 5.12e % 5.12e\n"%(atom[0],atom[1],atom[2])
        
    poscar=open(poscarName,"w")
    poscar.write(data)
    return poscar
