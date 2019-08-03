import math
import numpy as np
from model import *
spaceFrac=0.4
def gen3d(pbIBondLength,cation,metal,halide,alpha,beta,gamma,ligand):
    allAtoms=[]
    angle1=(alpha)*math.pi/180.0
    angle2=(beta)*math.pi/180.0
    angle3=(gamma)*math.pi/180.0


    rotx=np.array([[1.0,0.0,0.0],[0.0,math.cos(angle1),-1*math.sin(angle1)],[0.0,math.sin(angle1),math.cos(angle1)]])
    roty=np.array([[math.cos(angle2),0.0,math.sin(angle2)],[0.0,1.0,0.0],[-1*math.sin(angle2),0.0,math.cos(angle2)]])
    rotz=np.array([[math.cos(angle3),-1*math.sin(angle3),0.0],[math.sin(angle3),math.cos(angle3),0.0],[0.0,0.0,1.0]])
    pb1mat=np.array([[0.0,0.0,0.0]])
    i1mat=np.array([0.0,pbIBondLength,0.0])
    i2mat=np.array([pbIBondLength,0.0,0.0])
    i3mat=np.array([0.0,0.0,pbIBondLength])

    i1mat=i1mat.dot(rotx).dot(roty).dot(rotz)
    i2mat=i2mat.dot(rotx).dot(roty).dot(rotz)
    i3mat=i3mat.dot(rotx).dot(roty).dot(rotz)

    acell=4.0*i2mat[0]
    bcell=4.0*i1mat[1]
    ccell=4.0*i3mat[2]
    print('acell={:f}\n'.format(acell))
    print('bcell={:f}\n'.format(bcell))
    print('ccell={:f}\n'.format(ccell))

    acenter=0.25*acell
    bcenter=0.25*bcell
    ccenter=0.25*ccell
    modifiers=[[0.0, 0.0,0.0,1,1,1,'a'],\
        [0.0,0.0,0.5*ccell,1,1,-1,'b'],\
        [0.0,0.5*bcell,0.0,1,-1,1,'c'],\
        [0.0,0.5*bcell,0.5*ccell,1,-1,-1,'d'],\
        [0.5*acell, 0.0,0.0,-1,1,1,'e'],\
        [0.5*acell,0.0,0.5*ccell,-1,1,-1,'f'],\
        [0.5*acell,0.5*bcell,0.0,-1,-1,1,'g'],\
        [0.5*acell,0.5*bcell,0.5*ccell,-1,-1,-1,'h']]
    for m in modifiers:
        mode=m[3]*m[4]*m[5]
        allAtoms.append(Atom(metal+'1'+m[6],m[0],m[1],m[2],metal))
        allAtoms.append(Atom(halide+'1'+m[6],m[0]+mode*(i1mat[0]),m[1]+(i1mat[1]),m[2]+mode*(i1mat[2]),halide))
        allAtoms.append(Atom(halide+'2'+m[6],m[0]+(i2mat[0]),m[1]+mode*(i2mat[1]),m[2]+mode*(i2mat[2]),halide))
        allAtoms.append(Atom(halide+'3'+m[6],m[0]+mode*(i3mat[0]),m[1]+mode*(i3mat[1]),m[2]+(i3mat[2]),halide))
        for a in ligand.atomList:
            allAtoms.append(Atom(a.name+m[6],m[0]+acenter+a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))

    return Structure(name=cation+metal+halide,atomList=allAtoms,acell=acell,bcell=bcell,ccell=ccell)

def genRP(n,layerSpacing,pbIBondLength,cation,metal,halide,atilt,btilt,ctilt,ligand):
    allAtoms=[]
    #in angstroms
    print(layerSpacing)
    print(pbIBondLength)
    alpha=atilt
    beta=btilt
    gamma=ctilt
    #in degrees
    angle1=(alpha)*math.pi/180.0
    angle2=(beta)*math.pi/180.0
    angle3=(gamma)*math.pi/180.0

    rotx=np.array([[1.0,0.0,0.0],[0.0,math.cos(angle1),-1*math.sin(angle1)],[0.0,math.sin(angle1),math.cos(angle1)]])
    roty=np.array([[math.cos(angle2),0.0,math.sin(angle2)],[0.0,1.0,0.0],[-1*math.sin(angle2),0.0,math.cos(angle2)]])
    rotz=np.array([[math.cos(angle3),-1*math.sin(angle3),0.0],[math.sin(angle3),math.cos(angle3),0.0],[0.0,0.0,1.0]])
    pb1mat=np.array([[0.0,0.0,0.0]])
    i1mat=np.array([0.0,pbIBondLength,0.0])
    i2mat=np.array([pbIBondLength,0.0,0.0])
    i3mat=np.array([0.0,0.0,pbIBondLength])
    i4mat=np.array([-1*pbIBondLength,0.0,0.0])

    i1mat=i1mat.dot(rotx).dot(roty).dot(rotz)
    i2mat=i2mat.dot(rotx).dot(roty).dot(rotz)
    i3mat=i3mat.dot(rotx).dot(roty).dot(rotz)
    i4mat=i4mat.dot(rotx).dot(roty).dot(rotz)

    acell=layerSpacing
    bcell=4.0*i1mat[1]
    ccell=4.0*i3mat[2]
    print('acell={:f}\n'.format(acell))
    print('bcell={:f}\n'.format(bcell))
    print('ccell={:f}\n'.format(ccell))

    modifiers=[[0.0, 0.0,0.0,1,1,1,'a'],\
        [0.0,0.0,0.5*ccell,1,1,-1,'b'],\
        [0.0,0.5*bcell,0.0,1,-1,1,'c'],\
        [0.0,0.5*bcell,0.5*ccell,1,-1,-1,'d'],\
        [acell, 0.25*bcell,0.25*ccell,1,1,1,'e'],\
        [acell,0.25*bcell,0.75*ccell,1,1,-1,'f'],\
        [acell,0.75*bcell,0.25*ccell,1,-1,1,'g'],\
        [acell,0.75*bcell,0.75*ccell,1,-1,-1,'h']]
    bcenter=0.25*bcell
    ccenter=0.25*ccell
    print(bcenter)
    print(ccenter)
    for m in modifiers:
        mode=m[3]*m[4]*m[5]
        flip=-1
        for i in range(n):
            flip=flip*-1
            t=flip*mode
            allAtoms.append(Atom(metal+'1'+m[6],m[0]+2*i2mat[0]*i,m[1],m[2],metal))
            allAtoms.append(Atom(halide+'1'+m[6],m[0]+t*i1mat[0]+2*i2mat[0]*i,i1mat[1]+m[1],t*i1mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'2'+m[6],m[0]+i2mat[0]+2*i2mat[0]*i,t*i2mat[1]+m[1],t*i2mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'3'+m[6],m[0]+t*i3mat[0]+2*i2mat[0]*i,t*i3mat[1]+m[1],i3mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'4'+m[6],m[0]+i4mat[0]+2*i2mat[0]*i,t*i4mat[1]+m[1],t*i4mat[2]+m[2],halide))
        for j in range(n-1):
            acenter=(1+2*j)*i2mat[0]
            for a in cation.atomList:
                allAtoms.append(Atom(a.name+m[6],m[0]+acenter+a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))
        acenter=((n-1)*2)*i2mat[0]*(1-spaceFrac)+layerSpacing*spaceFrac
        for a in ligand.atomList:
            allAtoms.append(Atom(a.name+m[6],m[0]+acenter+a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))
        acenter=-spaceFrac*(layerSpacing-((n-1)*2)*i2mat[0])
        for a in ligand.atomList:
            allAtoms.append(Atom(a.name+m[6],m[0]+acenter-a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))
    return Structure(name='RP '+str(n)+metal+halide,atomList=allAtoms,acell=acell*2,bcell=bcell,ccell=ccell)

def genDJ(n,layerSpacing,pbIBondLength,cation,metal,halide,atilt,btilt,ctilt,ligand):
    allAtoms=[]
    #in angstroms
    print(layerSpacing)
    print(pbIBondLength)
    alpha=atilt
    beta=btilt
    gamma=ctilt
    #in degrees
    angle1=(alpha)*math.pi/180.0
    angle2=(beta)*math.pi/180.0
    angle3=(gamma)*math.pi/180.0

    rotx=np.array([[1.0,0.0,0.0],[0.0,math.cos(angle1),-1*math.sin(angle1)],[0.0,math.sin(angle1),math.cos(angle1)]])
    roty=np.array([[math.cos(angle2),0.0,math.sin(angle2)],[0.0,1.0,0.0],[-1*math.sin(angle2),0.0,math.cos(angle2)]])
    rotz=np.array([[math.cos(angle3),-1*math.sin(angle3),0.0],[math.sin(angle3),math.cos(angle3),0.0],[0.0,0.0,1.0]])
    pb1mat=np.array([[0.0,0.0,0.0]])
    i1mat=np.array([0.0,pbIBondLength,0.0])
    i2mat=np.array([pbIBondLength,0.0,0.0])
    i3mat=np.array([0.0,0.0,pbIBondLength])
    i4mat=np.array([-1*pbIBondLength,0.0,0.0])

    i1mat=i1mat.dot(rotx).dot(roty).dot(rotz)
    i2mat=i2mat.dot(rotx).dot(roty).dot(rotz)
    i3mat=i3mat.dot(rotx).dot(roty).dot(rotz)
    i4mat=i4mat.dot(rotx).dot(roty).dot(rotz)

    acell=layerSpacing
    bcell=4.0*i1mat[1]
    ccell=4.0*i3mat[2]
    print('acell={:f}\n'.format(acell))
    print('bcell={:f}\n'.format(bcell))
    print('ccell={:f}\n'.format(ccell))

    modifiers=[[0.0, 0.0,0.0,1,1,1,'a'],\
        [0.0,0.0,0.5*ccell,1,1,-1,'b'],\
        [0.0,0.5*bcell,0.0,1,-1,1,'c'],\
        [0.0,0.5*bcell,0.5*ccell,1,-1,-1,'d']]
    bcenter=0.25*bcell
    ccenter=0.25*ccell
    print(bcenter)
    print(ccenter)
    for m in modifiers:
        mode=m[3]*m[4]*m[5]
        flip=-1
        for i in range(n):
            flip=flip*-1
            t=flip*mode
            allAtoms.append(Atom(metal+'1'+m[6],m[0]+2*i2mat[0]*i,m[1],m[2],metal))
            allAtoms.append(Atom(halide+'1'+m[6],m[0]+t*i1mat[0]+2*i2mat[0]*i,i1mat[1]+m[1],t*i1mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'2'+m[6],m[0]+i2mat[0]+2*i2mat[0]*i,t*i2mat[1]+m[1],t*i2mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'3'+m[6],m[0]+t*i3mat[0]+2*i2mat[0]*i,t*i3mat[1]+m[1],i3mat[2]+m[2],halide))
            allAtoms.append(Atom(halide+'4'+m[6],m[0]+i4mat[0]+2*i2mat[0]*i,t*i4mat[1]+m[1],t*i4mat[2]+m[2],halide))
        for j in range(n-1):
            acenter=(2*j+1)*i2mat[0]
            for a in cation.atomList:
                allAtoms.append(Atom(a.name+m[6],m[0]+acenter+a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))
        acenter=(((n-1)*2)*i2mat[0]+layerSpacing)/2.0
        for a in ligand.atomList:
            allAtoms.append(Atom(a.name+m[6],m[0]+acenter+a.x,bcenter+m[1]+a.y,ccenter+m[2]+a.z,a.element))
    return Structure(name='DJ'+str(n)+metal+halide,atomList=allAtoms,acell=acell,bcell=bcell,ccell=ccell)
