import numpy as np
from tqdm import tqdm
from pymech.exadata import exadata, elem
import pandas as pd
import quadpy
from pymech.neksuite import readnek
from pymech.neksuite import writenek

def stacked_sum(var, (z,y,x)):
        return np.sum(np.sum(np.sum(var*z,axis=0)*y,axis=0)*x)

filename='..\Mesh0110_1p50.f00008'
field=readnek(filename)

nel=field.nel

lr1=field.lr1

scheme = quadpy.c1.gauss_lobatto(field.lr1[2])
gausslobattopoints=scheme.points
gausslobattoweights=scheme.weights
w=gausslobattoweights

TotalVol=0.0
TotalMass=0.0
TotalPotentialEnergy=0.0
TotalKineticEnergy=0.0

for n in np.arange(nel):
        ax=field.elem[n].pos[0,0,0,0]
        bx=field.elem[n].pos[0,0,0,-1]
        x=w*((bx-ax)/2.0)

        ay=field.elem[n].pos[1,0,0,0]
        by=field.elem[n].pos[1,0,-1,0]
        wyprime=w*((by-ay)/2.0)    

        az=field.elem[n].pos[2,0,0,0]
        bz=field.elem[n].pos[2,-1,0,0]
        wzprime=w*((bz-az)/2.0)

        temp=field.elem[n].vel[0,:,:,:]
        kineticenergy = 0.5*sum([field.elem[n].vel[i,:,:,:]**2 for i in [0,1,2]])

        zyx = (wzprime[:,np.newaxis,np.newaxis], wyprime[:,np.newaxis], wxprime)

        TotalVol += stacked_sum(np.ones(temp.shape), zyx)
        TotalMass += stacked_sum(field.elem[n].temp[0,:,:,:], zyx)       
        TotalPotentialEnergy += stacked_sum(field.elem[n].temp[0,:,:,:]*field.elem[n].pos[1,:,:,:], zyx)        
        TotalKineticEnergy += stacked_sum(kineticenergy, zyx)

 

print(TotalVol)
print(TotalMass)
print(TotalPotentialEnergy)
print(TotalKineticEnergy)