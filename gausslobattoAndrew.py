#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:13:49 2019
@author: andrewooi
"""
import os
import numpy as np
#from exadata import exadata
import pandas as p
import quadpy
from neksuite import readnek
from neksuite import writenek

os.chdir('D:\\Data\\aooi\\Nek5000\\StratifiedGravityCurrent\\NewMeshes\\Mesh0110_1p5\\')
filename='Mesh0110_1p50.f00008'
nzplanes=10 # number of macro z planes
data=readnek(filename)
nel=data.nel
xdata=np.zeros(nel)
ydata=np.zeros(nel)
zdata=np.zeros(nel)
index=np.zeros(nel,dtype=int)
lr1=data.lr1
zaveragedu=np.zeros((lr1[0],lr1[1]))
zaveragedv=np.zeros((lr1[0],lr1[1]))
zaveragedw=np.zeros((lr1[0],lr1[1]))
zaveragedpres=np.zeros((lr1[0],lr1[1]))
zaveragedtemp=np.zeros((lr1[0],lr1[1]))

scheme = quadpy.c1.gauss_lobatto(data.lr1[2])
gausslobattopoints=scheme.points
gausslobattoweights=scheme.weights

xi=gausslobattopoints
w=gausslobattoweights
#zaveragedu=np.zeros((lr1[1],lr1[0]))
element_number=np.arange(0,dtype=np.int32)
le=0

 

#

#Finding centroid of the macro elements to sort the data

#

for n in np.arange(nel):
    index[n]=int(n)
    xdata[n]=np.around((data.elem[n].pos[0,0,0,0]+data.elem[n].pos[0,0,0,lr1[0]-1])/2.0,decimals=3)
    ydata[n]=np.around((data.elem[n].pos[1,0,0,0]+data.elem[n].pos[1,0,lr1[0]-1,0])/2.0,decimals=3)
    zdata[n]=np.around((data.elem[n].pos[2,0,0,0]+data.elem[n].pos[2,lr1[0]-1,0,0])/2.0,decimals=3)

# Put data in pandas dataframe then sort the data by x, then y then z.
dataset = pd.DataFrame({'index': index[:],'xval': xdata[:],'yval': ydata[:], 'zval': zdata[:]})
sorteddataset=dataset.sort_values(by=['xval','yval','zval'])
#sorteddataset=dataset.sort_values(by=['yval','xval','zval'])

nplanarelements=int(nel/nzplanes)

for m in np.arange(nplanarelements):
#for m in np.arange(5):
    print(m)
    zaveragedu[:,:]=0.0
    zaveragedv[:,:]=0.0
    zaveragedw[:,:]=0.0
    zaveragedpres[:,:]=0.0
    zaveragedtemp[:,:]=0.0
    totallengthofz=0.0

    for n in np.arange(nzplanes):
        i=sorteddataset['index'].iloc[m*nzplanes+n]
        a=data.elem[i].pos[2,0,0,0]
        b=data.elem[i].pos[2,data.lr1[2]-1,0,0]
        z=((b-a)/2.0)*xi+(a+b)/2.0
        wprime=w*((b-a)/2.0)

        tempu=np.sum(data.elem[i].vel[0,:,:,:]*wprime[:,np.newaxis,np.newaxis],axis=0)
        totallengthofz=totallengthofz+np.sum(wprime)
        zaveragedu=zaveragedu+tempu

        tempv=np.sum(data.elem[i].vel[1,:,:,:]*wprime[:,np.newaxis,np.newaxis],axis=0)
        zaveragedv=zaveragedv+tempv

        tempw=np.sum(data.elem[i].vel[2,:,:,:]*wprime[:,np.newaxis,np.newaxis],axis=0)
        zaveragedw=zaveragedw+tempw

        temppres=np.sum(data.elem[i].pres[0,:,:,:]*wprime[:,np.newaxis,np.newaxis],axis=0)
        zaveragedpres=zaveragedpres+temppres

        temptemp=np.sum(data.elem[i].temp[0,:,:,:]*wprime[:,np.newaxis,np.newaxis],axis=0)
        zaveragedtemp=zaveragedtemp+temptemp

    zaveragedu=zaveragedu/totallengthofz
    zaveragedv=zaveragedv/totallengthofz
    zaveragedw=zaveragedw/totallengthofz
    zaveragedpres=zaveragedpres/totallengthofz
    zaveragedtemp=zaveragedtemp/totallengthofz

 

    # for n in np.arange(nzplanes):

    #     i=sorteddataset['index'].iloc[m*nzplanes+n]

    #     for p in np.arange(lr1[2]):

    #         data.elem[i].vel[0,p,:,:]=zaveragedu

    #         data.elem[i].vel[1,p,:,:]=zaveragedv

    #         data.elem[i].vel[2,p,:,:]=zaveragedw

    #         data.elem[i].pres[0,:,:,:]=zaveragedpres

    #         data.elem[i].temp[0,:,:,:]=zaveragedtemp

    for n in np.arange(nzplanes):
        i=sorteddataset['index'].iloc[m*nzplanes+n]
        data.elem[i].vel[0,:,:,:]=zaveragedu
        data.elem[i].vel[1,:,:,:]=zaveragedv
        data.elem[i].vel[2,:,:,:]=zaveragedw
        data.elem[i].pres[0,:,:,:]=zaveragedpres
        data.elem[i].temp[0,:,:,:]=zaveragedtemp

writenek('spanwiseavgMesh0110_1p50.f00008',data)