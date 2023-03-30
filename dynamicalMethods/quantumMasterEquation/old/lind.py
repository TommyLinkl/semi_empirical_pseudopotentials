import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA
import lind



def Lind_1lead(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
		        d1up, d1down, d2up,d2down,Ham, Rho,GammaL1,GammaL2): 

        Rho011 =  GammaL1*np.einsum('ba,bc,cd->ad',d1up,D1upLtilde,Rho, optimize = True)
        Rho011 += GammaL1* np.einsum('ba,bc,cd->ad',d1down,D1downLtilde,Rho, optimize = True)
        Rho011 += GammaL2* np.einsum('ba,bc,cd->ad',d2up,D2upLtilde,Rho, optimize = True)
        Rho011 += GammaL2*np.einsum('ba,bc,cd->ad',d2down,D2downLtilde,Rho, optimize = True)

        Rho011 +=  GammaL1* np.einsum('ab,cb,cd->ad',d1up,D1upL,Rho, optimize = True)
        Rho011 +=  GammaL1*np.einsum('ab,cb,cd->ad',d1down,D1downL,Rho, optimize = True)
        Rho011 +=  GammaL2*np.einsum('ab,cb,cd->ad',d2up,D2upL,Rho, optimize = True)
        Rho011 +=  GammaL2*np.einsum('ab,cb,cd->ad',d2down,D2downL,Rho, optimize = True)

        Rho011 -=  GammaL1*np.einsum('ba,bc,cd->ad',d1up,Rho,D1upL, optimize = True)
        Rho011 -=  GammaL1*np.einsum('ba,bc,cd->ad',d1down,Rho,D1downL, optimize = True)
        Rho011 -=  GammaL2*np.einsum('ba,bc,cd->ad',d2up,Rho,D2upL, optimize = True)
        Rho011 -=  GammaL2*np.einsum('ba,bc,cd->ad',d2down,Rho,D2downL, optimize = True)

        Rho011 -=  GammaL1*np.einsum('ab,bc,dc->ad',d1up,Rho,D1upLtilde, optimize = True)
        Rho011 -=  GammaL1*np.einsum('ab,bc,dc->ad',d1down,Rho,D1downLtilde, optimize = True)
        Rho011 -=  GammaL2*np.einsum('ab,bc,dc->ad',d2up,Rho,D2upLtilde, optimize = True)
        Rho011 -=  GammaL2*np.einsum('ab,bc,dc->ad',d2down,Rho,D2downLtilde, optimize = True)

        Rho011 = (Rho011 + np.conj(np.transpose(Rho011)))/2.0

        return -Rho011


def RK4_init(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
        	D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
        	d1up, d1down, d2up,d2down,Ham, Rho, GammaL1,GammaL2,GammaR1,GammaR2): 
    
        Rho0 = -1j* np.einsum('ab,bc->ac',Ham, Rho, optimize = True)
        Rho0 = (Rho0 + np.conj(np.transpose(Rho0)))/2.0

        Rho011L = Lind_1lead(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
        d1up, d1down, d2up,d2down,Ham, Rho, GammaL1,GammaL2)

        Rho011R = Lind_1lead(D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
        d1up, d1down, d2up,d2down,Ham, Rho, GammaR1,GammaR2)
        
        return Rho0 + Rho011L+ Rho011R

def RK4(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
        D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
        d1up, d1down, d2up,d2down,Ham, Rho, dt, GammaL1,GammaL2,GammaR1,GammaR2): 

        Rho_01 = lind.RK4_init(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
                D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
                d1up, d1down, d2up,d2down,Ham, Rho, GammaL1,GammaL2,GammaR1,GammaR2) 
        
        Rho_02 = lind.RK4_init(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
                D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
                d1up, d1down, d2up,d2down,Ham, Rho+Rho_01*dt/2.0,GammaL1,GammaL2,GammaR1,GammaR2) 
        
        Rho_03 = lind.RK4_init(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
                D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
                d1up, d1down, d2up,d2down,Ham, Rho+Rho_02*dt/2.0,GammaL1,GammaL2,GammaR1,GammaR2) 
        
        Rho_04 = lind.RK4_init(D1upLtilde,D1downLtilde,D2upLtilde,D2downLtilde,D1upL,D1downL,D2upL,D2downL,\
                D1upRtilde,D1downRtilde,D2upRtilde,D2downRtilde,D1upR,D1downR,D2upR,D2downR,\
                d1up, d1down, d2up,d2down,Ham, Rho+Rho_03*dt, GammaL1,GammaL2,GammaR1,GammaR2) 
        
        Rho = Rho + dt*(Rho_01 + 2.0*Rho_02 + 2.0*Rho_03 + Rho_04)/6.0

        return Rho
        

