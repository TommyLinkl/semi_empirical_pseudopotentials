###################################################################################
# External libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA

###################################################################################
# Internal files

import lind

###################################################################################
# Functions

# replace L by R for the other lead, generic for any lead 
def Lind_1lead(D_level_Ltilde, D_level_L, d_level, Ham, Rho, GammaL):
        
    dim = Rho.shape[0] 
    level = d_level.shape[0] 
    spin = d_level.shape[1] 
    
    Rho011 = np.zeros((dim, dim), dtype='complex128')
    
    for i in range(level):
        for sigma in range(spin):
            Rho011 +=  GammaL[i]*np.einsum('ba,bc,cd->ad', d_level[i, sigma], D_level_Ltilde[i, sigma], Rho, optimize=True)
            Rho011 +=  GammaL[i]*np.einsum('ab,cb,cd->ad', d_level[i, sigma], D_level_L[i, sigma], Rho, optimize=True)
            Rho011 -=  GammaL[i]*np.einsum('ba,bc,cd->ad', d_level[i, sigma], Rho, D_level_L[i, sigma], optimize=True)
            Rho011 -=  GammaL[i]*np.einsum('ab,bc,dc->ad', d_level[i, sigma], Rho, D_level_Ltilde[i, sigma], optimize=True)

    Rho011 = (Rho011 + np.conj(np.transpose(Rho011)))/2.0

    return -Rho011

###################################################################################
#

def RK4_init(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho, GammaL, GammaR): 

    Rho0 = -1j* np.einsum('ab,bc->ac', Ham, Rho, optimize=True)
    
    Rho0 = (Rho0 + np.conj(np.transpose(Rho0)))/2.0

    Rho011L = lind.Lind_1lead(D_level_Ltilde, D_level_L, d_level, Ham, Rho, GammaL)
    
    Rho011R = lind.Lind_1lead(D_level_Rtilde, D_level_R, d_level, Ham, Rho, GammaR)

    return Rho0 + Rho011L+ Rho011R

###################################################################################
#

def RK4(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho, dt, GammaL, GammaR): 

    Rho_01 = lind.RK4_init(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho, GammaL, GammaR)

    Rho_02 = lind.RK4_init(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho+Rho_01*dt/2.0, GammaL, GammaR)

    Rho_03 = lind.RK4_init(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho+Rho_02*dt/2.0, GammaL, GammaR)

    Rho_04 = lind.RK4_init(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho+Rho_03*dt, GammaL, GammaR)
                
    Rho = Rho + dt*(Rho_01 + 2.0*Rho_02 + 2.0*Rho_03 + Rho_04)/6.0

    return Rho
        
###################################################################################