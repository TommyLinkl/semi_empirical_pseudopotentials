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
#

def RK4_init(Ham, Rho): 

    Rho0 = -1j* np.einsum('ab,bc->ac', Ham, Rho )
    
    Rho0 = (Rho0 + np.conj(np.transpose(Rho0)))/2.0

    return Rho0 

###################################################################################
#

def RK4(Ham, Rho, dt): 

    Rho_01 = lind.RK4_init(Ham, Rho)

    Rho_02 = lind.RK4_init(Ham, Rho+Rho_01*dt/2.0)

    Rho_03 = lind.RK4_init(Ham, Rho+Rho_02*dt/2.0)

    Rho_04 = lind.RK4_init(Ham, Rho+Rho_03*dt)
                
    Rho = Rho + dt*(Rho_01 + 2.0*Rho_02 + 2.0*Rho_03 + Rho_04)/6.0

    return Rho
        
###################################################################################
