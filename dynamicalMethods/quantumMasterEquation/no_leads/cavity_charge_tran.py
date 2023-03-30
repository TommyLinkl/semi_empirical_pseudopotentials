import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA

import lind


level = 2
sites = 20
dim = sites

n_photon = 2


g01 = 0.01
hw = 0.09
kT=0.001
t = -0.00
dt = 1


E_level = np.zeros(level)
W_Coulomb = np.zeros((level,level,level,level))

E_level[0] = 0.01
E_level[1] = 0.1

Ham = np.zeros((dim,dim))


Ham[0,0] = E_level[0]+hw
for i in range(1,sites):
    Ham[i,i] = E_level[1]
    Ham[0,i] = g01 
    Ham[i,0] = g01 

for i in range(1,sites-1):
    Ham[i,i+1] = t 
    Ham[i+1,i] = t 

print(Ham)

evals, evecs = LA.eigh(Ham)
print(evals)
print(evecs)

Rho = np.zeros((dim,dim),dtype='complex128')
Rho[1,1] = 1.0 

itime=5400
y = np.zeros(itime)
for t in range(itime): 
    Rho = lind.RK4(Ham, Rho,dt)
    y[t] = Rho[-1,-1]
    
#plt.plot(bias,I_current)
#print(GammaL)
#print(D_level_Ltilde[0,1])
plt.plot(np.arange(itime),y)
plt.show()

    
    
