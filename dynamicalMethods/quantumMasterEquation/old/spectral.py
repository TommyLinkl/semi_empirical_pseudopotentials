import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA

import lind 


level = 2
dim = 4**level
mL = 0.04
mR = -0.04

dt = 1

E1 = 0.01
U1 = 0.03
GammaL1=0.01
GammaR1=0.01


E2 = 0.01
U2 = 0.00
U12 = 0.00
GammaL2=0.0
GammaR2=0.0


kT=0.001



Gamma1 = GammaL1+GammaR1
Gamma2 = GammaL2+GammaR2

Ham = np.zeros((dim,dim))

iden4 = np.identity(4)
#dup = np.zeros((4,4))
#ddown = np.zeros((4,4))

#dup[0,1] = 1.0 
#dup[2,3] = -1.0 

#ddown[0,2] = 1.0 
#ddown[1,3] = 1.0 

d1up = np.zeros((dim,dim))
d1down = np.zeros((dim,dim))

d2up = np.zeros((dim,dim))
d2down = np.zeros((dim,dim))

num_operator_single = np.zeros((4,4))
num_operator_total = np.zeros((dim,dim))

Rho = np.zeros((dim,dim),dtype='complex128')
Rho[0,0] = 1.0 

bias=np.arange(0.0,0.12,0.006)

def fermiL(x,mL):
    return 1.0/(np.exp((x-mL)/kT)+1.0) 

def fermiR(x,mR):
    return 1.0/(np.exp((x-mR)/kT)+1.0) 


num_operator_single[1,1]=1.0
num_operator_single[2,2]=1.0
num_operator_single[3,3]=2.0

num_operator_total = np.kron(num_operator_single, iden4)+ np.kron(iden4, num_operator_single) 

Ham[0,0] = 0.0 
Ham[1,1] = E1 
Ham[2,2] = E1 
Ham[3,3] = 2.0*E1+U1 

Ham[4,4] = E2 
Ham[5,5] = E1+E2 + U12 
Ham[6,6] = E1+E2 + U12 
Ham[7,7] = 2*E1+U1+E2 + 2.0*U12 

Ham[8,8] = E2 
Ham[9,9] = E1+E2 + U12 
Ham[10,10] = E1+E2 + U12 
Ham[11,11] = 2*E1+U1+E2 + 2.0*U12 

Ham[12,12] = 2*E2+U2 
Ham[13,13] = E1+2*E2+U2 +  2.0*U12 
Ham[14,14] = E1+2*E2+U2 + 2.0*U12 
Ham[15,15] = 2*E1+U1+2*E2+U2 + 4.0*U12 


for i in range(4**(level-1)):
    d1up[0+4*i,1+4*i] = 1.0 
    d1up[2+4*i,3+4*i] = -1.0 

    d1down[0+4*i,2+4*i] = 1.0 
    d1down[1+4*i,3+4*i] = 1.0 

for i in range(4):
    d2up[i,1*4+i] = 1.0 
    d2up[2*4+i,3*4+i] = -1.0 

    d2down[i,2*4+i] = 1.0 
    d2down[1*4+i,3*4+i] = 1.0 


evals, evecs = LA.eigh(Ham)

D1up = np.einsum('ba,bc,cd->ad', evecs,d1up,evecs, optimize=True)
D1down = np.einsum('ba,bc,cd->ad', evecs,d1down,evecs, optimize=True)
D2up = np.einsum('ba,bc,cd->ad', evecs,d2up,evecs, optimize=True)
D2down = np.einsum('ba,bc,cd->ad', evecs,d2down,evecs, optimize=True)

D1upLtilde= np.zeros((dim,dim))
D1downLtilde= np.zeros((dim,dim))
D2upLtilde= np.zeros((dim,dim))
D2downLtilde= np.zeros((dim,dim))

D1upL= np.zeros((dim,dim))
D1downL= np.zeros((dim,dim))
D2upL= np.zeros((dim,dim))
D2downL= np.zeros((dim,dim))

D1upRtilde= np.zeros((dim,dim))
D1downRtilde= np.zeros((dim,dim))
D2upRtilde= np.zeros((dim,dim))
D2downRtilde= np.zeros((dim,dim))

D1upR= np.zeros((dim,dim))
D1downR= np.zeros((dim,dim))
D2upR= np.zeros((dim,dim))
D2downR = np.zeros((dim,dim))

I_current = np.zeros(bias.shape[0])

for i_bias in range(bias.shape[0]): 
    mL = -bias[i_bias]/2.0
    mR = bias[i_bias]/2.0
    for i in range(dim):
        for j in range(dim):
            D1upLtilde[i,j] = D1up[i,j]*(1.0-fermiL(evals[j]-evals[i],mL))
            D1downLtilde[i,j] = D1down[i,j]*(1.0-fermiL(evals[j]-evals[i],mL))
            D2upLtilde[i,j] = D2up[i,j]*(1.0-fermiL(evals[j]-evals[i],mL))
            D2downLtilde[i,j] = D2down[i,j]*(1.0-fermiL(evals[j]-evals[i],mL))
    
            D1upL[i,j] = D1up[i,j]*fermiL(evals[j]-evals[i],mL)
            D1downL[i,j] = D1down[i,j]*fermiL(evals[j]-evals[i],mL)
            D2upL[i,j] = D2up[i,j]*fermiL(evals[j]-evals[i],mL)
            D2downL[i,j] = D2down[i,j]*fermiL(evals[j]-evals[i],mL)
    
            D1upRtilde[i,j] = D1up[i,j]*(1.0-fermiR(evals[j]-evals[i],mR))
            D1downRtilde[i,j] = D1down[i,j]*(1.0-fermiR(evals[j]-evals[i],mR))
            D2upRtilde[i,j] = D2up[i,j]*(1.0-fermiR(evals[j]-evals[i],mR))
            D2downRtilde[i,j] = D2down[i,j]*(1.0-fermiR(evals[j]-evals[i],mR))
    
            D1upR[i,j] = D1up[i,j]*fermiR(evals[j]-evals[i],mR)
            D1downR[i,j] = D1down[i,j]*fermiR(evals[j]-evals[i],mR)
            D2upR[i,j] = D2up[i,j]*fermiR(evals[j]-evals[i],mR)
            D2downR[i,j] = D2down[i,j]*fermiR(evals[j]-evals[i],mR)
    
    
    D1upLtilde = np.einsum('ab,bc,dc->ad', evecs,D1upLtilde,evecs, optimize=True)
    D1downLtilde = np.einsum('ab,bc,dc->ad', evecs,D1downLtilde,evecs, optimize=True)
    D2upLtilde = np.einsum('ab,bc,dc->ad', evecs,D2upLtilde,evecs, optimize=True)
    D2downLtilde = np.einsum('ab,bc,dc->ad', evecs,D2downLtilde,evecs, optimize=True)
    
    D1upL = np.einsum('ab,bc,dc->ad', evecs,D1upL,evecs, optimize=True)
    D1downL = np.einsum('ab,bc,dc->ad', evecs,D1downL,evecs, optimize=True)
    D2upL = np.einsum('ab,bc,dc->ad', evecs,D2upL,evecs, optimize=True)
    D2downL = np.einsum('ab,bc,dc->ad', evecs,D2downL,evecs, optimize=True)
    
    D1upRtilde = np.einsum('ab,bc,dc->ad', evecs,D1upRtilde,evecs, optimize=True)
    D1downRtilde = np.einsum('ab,bc,dc->ad', evecs,D1downRtilde,evecs, optimize=True)
    D2upRtilde = np.einsum('ab,bc,dc->ad', evecs,D2upRtilde,evecs, optimize=True)
    D2downRtilde = np.einsum('ab,bc,dc->ad', evecs,D2downRtilde,evecs, optimize=True)
    
    D1upR = np.einsum('ab,bc,dc->ad', evecs,D1upR,evecs, optimize=True)
    D1downR = np.einsum('ab,bc,dc->ad', evecs,D1downR,evecs, optimize=True)
    D2upR = np.einsum('ab,bc,dc->ad', evecs,D2upR,evecs, optimize=True)
    D2downR = np.einsum('ab,bc,dc->ad', evecs,D2downR,evecs, optimize=True)
    
    
    pop1 = np.zeros(1000)
    pop2 = np.zeros(1000)
    
    for t in range(1000): 
        pop1[t] = np.real(Rho[1,1]+Rho[2,2]+2.0*Rho[3,3]+Rho[5,5]+Rho[6,6]+2.0*Rho[7,7]+Rho[9,9]+Rho[10,10]+2.0*Rho[11,11]+Rho[13,13]+Rho[14,14]+2.0*Rho[15,15])
        pop2[t] = np.real(Rho[4,4]+Rho[5,5]+Rho[6,6]+Rho[7,7]+Rho[8,8]+Rho[9,9]+Rho[10,10]+Rho[11,11]+2.0*Rho[12,12]+2.0*Rho[13,13]+2.0*Rho[14,14]+2.0*Rho[15,15])
        Rho = lind.RK4(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
            D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
            d1up, d1down, d2up,d2down,Ham, Rho, dt,GammaL1,GammaL2,GammaR1,GammaR2)  
    
    # include i[H, \rho] later
    I_current_max_1 = lind.Lind_1lead(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
        d1up, d1down, d2up,d2down,Ham, Rho, GammaL1,GammaL2)

    num_rho = np.einsum('ab,bc->ac', num_operator_total, Rho,optimize=True)

    # include i[H, num_rho] later
    I_current_max_2 = lind.Lind_1lead(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
        d1up, d1down, d2up,d2down,Ham, num_rho, GammaL1,GammaL2)

    I_current_1 = np.real(np.einsum('ab,ba->', num_operator_total, I_current_max_1, optimize=True))
    I_current_2 = np.real(np.einsum('aa->', I_current_max_2, optimize=True))

    I_current[i_bias] = I_current_1 - I_current_2



plt.plot(bias,I_current)
plt.show()

    
    
