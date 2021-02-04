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
        wxprime=w*((bx-ax)/2.0)

        ay=field.elem[n].pos[1,0,0,0]
        by=field.elem[n].pos[1,0,-1,0]
        wyprime=w*((by-ay)/2.0)    

        az=field.elem[n].pos[2,0,0,0]
        bz=field.elem[n].pos[2,-1,0,0]
        wzprime=w*((bz-az)/2.0)

        temp=field.elem[n].vel[0,:,:,:]
        kineticenergy=0.5*(field.elem[n].vel[0,:,:,:]*field.elem[n].vel[0,:,:,:]+field.elem[n].vel[1,:,:,:]*field.elem[n].vel[1,:,:,:]+field.elem[n].vel[2,:,:,:]*field.elem[n].vel[2,:,:,:])

        #tempu=np.sum(np.ones(temp.shape)*wzprime[:,np.newaxis,np.newaxis],axis=0)
        TotalVol=TotalVol+np.sum(np.sum(np.sum(np.ones(temp.shape)*wzprime[:,np.newaxis,np.newaxis],axis=0)*wyprime[:,np.newaxis],axis=0)*wxprime)     
        TotalMass=TotalMass+np.sum(np.sum(np.sum(field.elem[n].temp[0,:,:,:]*wzprime[:,np.newaxis,np.newaxis],axis=0)*wyprime[:,np.newaxis],axis=0)*wxprime)
        TotalPotentialEnergy=TotalPotentialEnergy+np.sum(np.sum(np.sum(field.elem[n].temp[0,:,:,:]*field.elem[n].pos[1,:,:,:]*wzprime[:,np.newaxis,np.newaxis],axis=0)*wyprime[:,np.newaxis],axis=0)*wxprime)
        TotalKineticEnergy=TotalKineticEnergy+np.sum(np.sum(np.sum(kineticenergy*wzprime[:,np.newaxis,np.newaxis],axis=0)*wyprime[:,np.newaxis],axis=0)*wxprime)

 

print(TotalVol)
print(TotalMass)
print(TotalPotentialEnergy)
print(TotalKineticEnergy)