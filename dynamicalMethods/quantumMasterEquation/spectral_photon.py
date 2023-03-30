import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA
import sys

import lind_2 


spin  = 2
level = 2
dim   = 4**level

n_photon = 10

dim_total = dim*n_photon 

g00 = 0.04
g01 = 0.00
hw  = 0.05
kT  = 0.004
dt  = 1
Ed  = ( g00*g00 / (hw) )

E_level = np.zeros(level)
W_Coulomb = np.zeros((level,level,level,level))

E_level[0] = Ed
E_level[1] = 10.0
W_Coulomb[0,0,0,0] = 0.00 

GammaL = np.zeros(level)
GammaR = np.zeros(level)

GammaL[0] = 0.01
GammaR[0] = 0.01
n_photonGamma = GammaL+GammaR

Ham = np.zeros((dim,dim))
Ham_total = np.zeros((dim_total,dim_total))
d_level = np.zeros((level,spin,dim,dim))

iden4 = np.identity(4)
iden4back = np.identity(4)
iden4back[1,1] = -1.0 
iden4back[2,2] = -1.0 
dup = np.zeros((4,4))
ddown = np.zeros((4,4))

dup[0,1] = 1.0 
dup[2,3] = 1.0 
ddown[0,2] = 1.0 
ddown[1,3] = -1.0 


for i in range(level):
   dumup = dup.copy()
   dumdown = ddown.copy()
   for j in range(i):
       dumup = np.kron(dumup,iden4)
       dumdown = np.kron(dumdown,iden4)

   for j in range(i+1,level):
       dumup = np.kron(iden4back, dumup)
       dumdown = np.kron(iden4back, dumdown)

   d_level[i,0] = dumup.copy()
   d_level[i,1] = dumdown.copy()

bias = np.arange(0.0, 0.201, 0.01)

def fermi(x,mu):
    return 1.0/(np.exp((x-mu)/kT)+1.0) 

num_operator_single = np.zeros((4,4))
num_operator_total = np.zeros((dim,dim))
num_operator_single[1,1]=1.0
num_operator_single[2,2]=1.0
num_operator_single[3,3]=2.0


for i in range(level):
   dum = num_operator_single.copy()
   for j in range(i):
       dum = np.kron(iden4,dum)
   for j in range(i+1,level):
       dum = np.kron(dum,iden4)
   num_operator_total += dum

spin_operator_single = np.zeros((4,4))
spin_operator_total = np.zeros((dim,dim))
spin_operator_single[1,1]=1.0
spin_operator_single[2,2]=-1.0
spin_operator_single[3,3]=0.0


for i in range(level):
   dum = spin_operator_single.copy()
   for j in range(i):
       dum = np.kron(iden4,dum)
   for j in range(i+1,level):
       dum = np.kron(dum,iden4)
   spin_operator_total += dum


for i in range(level):
    for sigma in range(spin):
        Ham += E_level[i]*np.einsum('ab,bc->ac',np.transpose(d_level[i,sigma]),d_level[i,sigma])

for i in range(level):
    for j in range(level):
        for k in range(level):
            for l in range(level):
                dum = np.einsum('ab,bc->ac',d_level[l,1],d_level[j,0])
                dum = np.einsum('ab,bc->ac',np.transpose(d_level[k,1]),dum)
                dum = np.einsum('ab,bc->ac',np.transpose(d_level[i,0]),dum)
                Ham += W_Coulomb[i,j,k,l]*dum


# define photon creation operator
a_creation = np.zeros((n_photon,n_photon))
iden_photon = np.identity(n_photon)
iden_electron = np.identity(dim)

d_0_up_total = np.kron(d_level[0,0],iden_photon)
d_0_down_total = np.kron(d_level[0,1],iden_photon)

d_1_up_total = np.kron(d_level[1,0],iden_photon)
d_1_down_total = np.kron(d_level[1,1],iden_photon)

d_level_total = np.zeros((level, spin, dim_total, dim_total))
for i in range(level):
    for sigma in range(spin):
        d_level_total[i, sigma] = np.kron(d_level[i, sigma], iden_photon)

for a in range(n_photon-1):
    a_creation[a+1,a] = np.sqrt(1.0*(a+1.0)) 
    
a_creation_total = np.kron(iden_electron,a_creation)

x_total = a_creation_total + np.transpose(a_creation_total)

for j in range(dim):
    for a in range(n_photon):
        Ham_total[j*n_photon+a,j*n_photon+a]= Ham[j,j] + hw*(a)#+0.5)

#dum = np.einsum('ab,bc->ac',d_0_up_total,np.transpose(d_1_up_total))
#dum += np.einsum('ab,bc->ac',d_0_down_total,np.transpose(d_1_down_total))
#dum += np.transpose(dum)
#Ham_total += g01*np.einsum('ab,bc->ac',dum, x_total)

dum = np.einsum('ab,bc->ac',np.transpose(d_0_up_total),d_0_up_total)
dum += np.einsum('ab,bc->ac',np.transpose(d_0_down_total),d_0_down_total)
#dum += np.transpose(dum)
Ham_total += g00*np.einsum('ab,bc->ac', dum, x_total)


np.savetxt("Ham.dat", Ham_total, fmt="% .3f", delimiter=" ")
np.savetxt("a_creation.dat", a_creation, fmt="% .3f", delimiter=" ")
np.savetxt("a_creation_total.dat", a_creation_total, fmt="% .3f", delimiter=" ")
np.savetxt("d0d1d1d0.dat", dum, fmt="% .3f", delimiter=" ")
np.savetxt("xOperator.dat", x_total, fmt="% .3f", delimiter=" ")
#sys.exit()

Ham = Ham_total
dim = dim_total
d_level = d_level_total

num_operator_total = np.kron(num_operator_total, iden_photon)

#print(d_level[0,0])
#print(d_level[0,1])
Rho = np.zeros((dim,dim),dtype='complex128')
Rho[0,0] = 1.0 

evals, evecs = LA.eigh(Ham)

np.savetxt("eigenenergies.dat", evals, fmt="% .6f", delimiter=" ")
np.savetxt("eigenstates.dat", evecs, fmt="% .3f", delimiter=" ")

D_level= np.zeros((level,spin,dim,dim))
for i in range(level):
    for sigma in range(spin):
        D_level[i,sigma] = np.einsum('ba,bc,cd->ad', evecs,d_level[i,sigma],evecs, optimize=True)

D_level_L= np.zeros((level,spin,dim,dim))
D_level_Ltilde= np.zeros((level,spin,dim,dim))
D_level_R= np.zeros((level,spin,dim,dim))
D_level_Rtilde= np.zeros((level,spin,dim,dim))

I_current = np.zeros(bias.shape[0])

for i_bias in range(bias.shape[0]): 
    mL = -bias[i_bias]/2.0
    mR = bias[i_bias]/2.0
    for i_level in range(level): 
        for sigma in range(spin): 
            for i in range(dim):
                for j in range(dim):
                    D_level_Ltilde[i_level,sigma,i,j] = D_level[i_level,sigma,i,j]*(1.0-fermi(evals[j]-evals[i],mL))
                    D_level_L[i_level,sigma,i,j] = D_level[i_level,sigma,i,j]*fermi(evals[j]-evals[i],mL)
    
                    D_level_Rtilde[i_level,sigma,i,j] = D_level[i_level,sigma,i,j]*(1.0-fermi(evals[j]-evals[i],mR))
                    D_level_R[i_level,sigma,i,j] = D_level[i_level,sigma,i,j]*fermi(evals[j]-evals[i],mR)
    
#    print(evecs.shape) 
    for i_level in range(level): 
        for sigma in range(spin): 
            D_level_Ltilde[i_level,sigma] = np.einsum('ab,bc,dc->ad', evecs,D_level_Ltilde[i_level,sigma],evecs, optimize=True)
            D_level_L[i_level,sigma]= np.einsum('ab,bc,dc->ad', evecs,D_level_L[i_level,sigma],evecs, optimize=True)
            D_level_R[i_level,sigma]= np.einsum('ab,bc,dc->ad', evecs,D_level_R[i_level,sigma],evecs, optimize=True)
            D_level_Rtilde[i_level,sigma] = np.einsum('ab,bc,dc->ad', evecs,D_level_Rtilde[i_level,sigma],evecs, optimize=True)
    
    
    pop_total = np.zeros(2000)
    spin_total = np.zeros(2000)
    rhoTrace = np.zeros(2000)

    for t in range(2000): 
        Rho = lind_2.RK4(D_level_Ltilde,D_level_L,D_level_Rtilde,D_level_R, d_level,Ham, Rho, dt,GammaL,GammaR)
        #spin_total[t] = np.einsum('ab,ba', spin_operator_total, Rho,optimize=True)
        if i_bias == (bias.shape[0]-1):
            rhoTrace[t] = np.einsum('ab,ba->', Rho, num_operator_total)
    
    if i_bias == (bias.shape[0]-1):
        np.savetxt("rhoTrace.dat", rhoTrace, fmt="% .6f", delimiter=" ")

    I_current_max_1 =  lind_2.Lind_1lead(D_level_Ltilde,D_level_L,d_level, Ham, Rho, GammaL)

    num_rho = np.einsum('ab,bc->ac', num_operator_total, Rho,optimize=True)

    I_current_max_2 =  lind_2.Lind_1lead(D_level_Ltilde,D_level_L,d_level, Ham,num_rho, GammaL)

    I_current_1 = np.real(np.einsum('ab,ba->', num_operator_total, I_current_max_1, optimize=True))
    I_current_2 = np.real(np.einsum('aa->', I_current_max_2, optimize=True))
    
    #plt.plot(np.arange(2000),spin_total)

    I_current[i_bias] = I_current_1 - I_current_2



#plt.plot(bias,I_current)
#print(GammaL)
#print(D_level_Ltilde[0,1])
#plt.plot(np.arange(2000),I_current)
#plt.show()
np.savetxt("I-V-curve.dat", I_current, fmt="% .6f", delimiter=" ")
#np.savetxt("pop-v-curve.dat", pop_total, fmt="% .6f", delimiter=" ")
#np.savetxt("spin-v-curve.dat", spin_total, fmt="% .6f", delimiter=" ")
    
    
