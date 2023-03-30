###################################################################################
"""
  
  Program that ...

  Authors: 
  Last Modified: 

  Input files:
    energies.par -> 
    coulomb.par ->
    gamma.par ->
    bias.par ->
    input.par -> 

  Output files:
    xx -> xx

  Example run commands:
    $ python spectral_python2.py --kT=0.01 --nTimeSteps=100 --lowestEnergyLevel=-0.2
    $ python spectral_python2.py @input.par

"""

###################################################################################
# External libraries

import sys
import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import linalg as LA
from argparse import ArgumentParser

###################################################################################
# Internal files

import lind_python2

###################################################################################
# Functions

def buildParser():
  """
  Parser for command line arguments
  Options:
    --outputDirNum=1
  """

  parser = ArgumentParser(fromfile_prefix_chars="@")

  parser.add_argument("--spinDegeneracy", type=int, default=2, 
                        help="the spin degeneracy of each single-particle state")
  
  parser.add_argument("--nTimeSteps", type=int, default=1000, help="the number of time steps")
  
  parser.add_argument("--dt", type=float, default=1.0, help="time step size")
  
  parser.add_argument("--kT", type=float, default=0.001, help="the temperature multiplied by Boltzmann constant")

  parser.add_argument("--lowestEnergyLevel", type=float, default=0.0, help="energy to shift the lowest single-particle state to")

  parser.add_argument("--biasEnergyCenter", type=float, default=0.0, help="energy of the middle of the bias")

  return parser

###################################################################################

def fermi(x, mu, kT):
  return 1.0/(np.exp((x-mu)/kT) + 1.0) 

###################################################################################

def readSingleParticleEnergyLevelInput(fileName):
  """
  Input:
    - 
  Returns:
    - spEnergyLevels = 
  """

  spEnergyLevels = np.loadtxt(fileName, dtype=float)

  return np.array(spEnergyLevels, ndmin=1)

###################################################################################

def readGammaInput(fileName):
  """
  Input:
    - 
  Returns:
    - gammaLeft =
    - gammaRight 
  """

  gamma = np.loadtxt(fileName, dtype=float, ndmin=2)

  gammaLeft  = np.copy(gamma[:, 0])
  gammaRight = np.copy(gamma[:, 1])

  return np.array(gammaLeft, ndmin=1), np.array(gammaRight, ndmin=1)

###################################################################################

def readBiasInput(fileName):
  """
  Input:
    - 
  Returns:
    - minBias =
    - maxBias
    - biasStepSize 
  """

  bias = np.loadtxt(fileName, dtype=float, ndmin=1)

  return bias[0], bias[1], bias[2]

###################################################################################

def readCoulombMatrixElements(nSingleParticleStates, fileName):
  """
  Input:
    - 
  Returns:
    - coulombMatrixElements =
  """

  n = nSingleParticleStates
  coulombMatrixElements = np.zeros((n, n, n, n), dtype=float)

  with open(fileName) as f:
    for line in f:
      row = line.split()
      coulombMatrixElements[int(row[0]), int(row[1]), int(row[2]), int(row[3])] = np.float_(row[4])

  return coulombMatrixElements

###################################################################################

def calcCurrent(D_level_Ltilde, D_level_L, d_level, num_operator_total, Ham, Rho, gammaLeft):
  """
  Input:
    - 
  Returns:
    - current =
  """
  
  I_current_max_1 = lind_python2.Lind_1lead(D_level_Ltilde, D_level_L, d_level, Ham, Rho, gammaLeft)

  num_rho = np.einsum('ab,bc->ac', num_operator_total, Rho)

  I_current_max_2 = lind_python2.Lind_1lead(D_level_Ltilde, D_level_L, d_level, Ham, num_rho, gammaLeft)

  I_current_1 = np.real(np.einsum('ab,ba->', num_operator_total, I_current_max_1))
  
  I_current_2 = np.real(np.einsum('aa->', I_current_max_2))

  return (I_current_1 - I_current_2)

###################################################################################

def main():
  """

  """
  ###################################################################################
  # Parameters

  parser = buildParser()
  params = parser.parse_args()

  ###################################################################################
  # Read in the input files used to build the Hamiltonian and Lindblad superoperators 

  # Single-particle energies 
  spEnergyLevels = readSingleParticleEnergyLevelInput("energies.par")
  energyLevelShift = (params.lowestEnergyLevel - spEnergyLevels[0])
  spEnergyLevels = spEnergyLevels + energyLevelShift
  np.savetxt("shiftedEnergies.dat", spEnergyLevels, fmt="% .12f", delimiter="\n")
  nLevels = spEnergyLevels.shape[0]
  dim = 4**nLevels

  # The Coulomb matrix elements between the single-particle states
  coulombMatrixElements = readCoulombMatrixElements(nLevels, "coulomb.par")

  # Coupling strengths between the electron baths each single-particle energy level
  gammaLeft, gammaRight = readGammaInput("gamma.par")
  gammaTotal = gammaLeft + gammaRight

  # Bias values to be studied
  minBias, maxBias, biasStepSize = readBiasInput("bias.par")
  bias = np.arange(minBias, maxBias, biasStepSize)

  ###################################################################################
  # Build the creation and annihilation operators for each level

  # Final annihilation operator will be stored here 
  d_level = np.zeros((nLevels, params.spinDegeneracy, dim, dim))

  # Identity operator for a single site
  iden4 = np.identity(4)
  iden4back = np.identity(4)
  iden4back[1, 1] = -1.0
  iden4back[2, 2] = -1.0

  # Creation and annihilation operator for a single site
  dup = np.zeros((4, 4))
  ddown = np.zeros((4, 4))
  dup[0, 1] = 1.0 
  dup[2, 3] = 1.0 
  ddown[0, 2] = 1.0 
  ddown[1, 3] = -1.0 

  # Build the total operator
  for i in range(nLevels):
    dumup = dup.copy()
    dumdown = ddown.copy()
    for j in range(i):
      dumup = np.kron(dumup, iden4)
      dumdown = np.kron(dumdown, iden4)
    for j in range(i+1, nLevels):
      dumup = np.kron(iden4back, dumup)
      dumdown = np.kron(iden4back, dumdown)
    d_level[i, 0] = dumup.copy()
    d_level[i, 1] = dumdown.copy()

  ###################################################################################
  # Build the total number operator

  num_operator_total = np.zeros((dim, dim))

  # Number operator for a single level
  num_operator_single = np.zeros((4, 4))
  num_operator_single[1, 1] = 1.0
  num_operator_single[2, 2] = 1.0
  num_operator_single[3, 3] = 2.0
  
  for i in range(nLevels):
    dum = num_operator_single.copy()
    for j in range(i):
      dum = np.kron(iden4, dum)
    for j in range(i+1, nLevels):
      dum = np.kron(dum, iden4)
    num_operator_total += dum

  ###################################################################################
  # Create total spin operators: Sx, Sy, Sz and S^2

  Sx_operator_total = np.zeros((dim, dim), dtype="complex128")
  Sy_operator_total = np.zeros((dim, dim), dtype="complex128")
  Sz_operator_total = np.zeros((dim, dim), dtype="complex128")

  # Spin operators for a single level
  Sx_operator_single = np.zeros((4, 4), dtype="complex128")
  Sx_operator_single[1, 2] = 1.0
  Sx_operator_single[2, 1] = 1.0

  Sy_operator_single = np.zeros((4, 4), dtype="complex128")
  Sy_operator_single[1, 2] = -1.0*1j
  Sy_operator_single[2, 1] =  1.0*1j

  Sz_operator_single = np.zeros((4, 4), dtype="complex128")
  Sz_operator_single[1, 1] =  1.0
  Sz_operator_single[2, 2] = -1.0
  
  # Build total Sx
  for i in range(nLevels):
    dum = Sx_operator_single.copy()
    for j in range(i):
      dum = np.kron(iden4, dum)
    for j in range(i+1, nLevels):
      dum = np.kron(dum, iden4)
    Sx_operator_total += dum

  # Build total Sy
  for i in range(nLevels):
    dum = Sy_operator_single.copy()
    for j in range(i):
      dum = np.kron(iden4, dum)
    for j in range(i+1, nLevels):
      dum = np.kron(dum, iden4)
    Sy_operator_total += dum

  # Build total Sz
  for i in range(nLevels):
    dum = Sz_operator_single.copy()
    for j in range(i):
      dum = np.kron(iden4, dum)
    for j in range(i+1, nLevels):
      dum = np.kron(dum, iden4)
    Sz_operator_total += dum

  # Build total S^2
  S2_operator_total = np.einsum('ab,bc->ac', Sx_operator_total, Sx_operator_total)
  S2_operator_total += np.einsum('ab,bc->ac', Sy_operator_total, Sy_operator_total)
  S2_operator_total += np.einsum('ab,bc->ac', Sz_operator_total, Sz_operator_total)

  ###################################################################################
  # Build the number operator and Hamiltonian

  # System Hamiltonian matrix
  Ham = np.zeros((dim, dim))

  for i in range(nLevels):
    for sigma in range(params.spinDegeneracy):
      Ham += spEnergyLevels[i]*np.einsum('ab,bc->ac', np.transpose(d_level[i, sigma]), d_level[i, sigma])

  for i in range(nLevels):
    for i_sigma in range(params.spinDegeneracy):
      for j in range(nLevels): # j_sigma == i_sigma
        for k in range(nLevels): 
          for k_sigma in range(params.spinDegeneracy):
            for l in range(nLevels): # l_sigma == k_sigma
              dum = np.einsum('ab,bc->ac', d_level[l, k_sigma], d_level[j, i_sigma])
              dum = np.einsum('ab,bc->ac', np.transpose(d_level[k, k_sigma]), dum)
              dum = np.einsum('ab,bc->ac', np.transpose(d_level[i, i_sigma]), dum)
              Ham += coulombMatrixElements[i, k, j, l]*dum/2.0

  # Print the initial system Hamiltonian
  np.savetxt("initialHamiltonian.dat", Ham, fmt="% .12f", delimiter=" ")

  ###################################################################################
  # Diagonalize the system Hamiltonian and store the eigenvalues and eigenvectors

  # Diagonalize the Hamiltonian using numpy's linear algebra library 
  evals, evecs = LA.eigh(Ham)

  # Print the system eigenvalues and eigenvectors
  np.savetxt("systemEigenvalues.dat", evals, fmt="% .12f", delimiter=" ")
  np.savetxt("systemEigenstates.dat", evecs, fmt="% .12f", delimiter=" ")

  ###################################################################################
  # Calculate <N>, <S^2> and <Sz> for each eigenstates of the system Hamiltonian
  
  S2_operator_eigenbasis = np.real(np.einsum('ba,bc,cd->ad', evecs, S2_operator_total, evecs))
  Sz_operator_eigenbasis = np.real(np.einsum('ba,bc,cd->ad', evecs, Sz_operator_total, evecs))
  N_operator_eigenbasis = np.real(np.einsum('ba,bc,cd->ad', evecs, num_operator_total, evecs))
  np.savetxt("S2_operator_eigenbasis.dat", S2_operator_eigenbasis, fmt="% .3f", delimiter=" ")
  np.savetxt("Sz_operator_eigenbasis.dat", Sz_operator_eigenbasis, fmt="% .3f", delimiter=" ")
  np.savetxt("N_operator_eigenbasis.dat", N_operator_eigenbasis, fmt="% .3f", delimiter=" ")

  # Print all quantum numbers for the system eigenstates
  np.savetxt("systemQuantumNumbers.dat", np.vstack((evals, np.diagonal(N_operator_eigenbasis), 0.25*np.diagonal(S2_operator_eigenbasis), 
	  												   0.5*np.diagonal(Sz_operator_eigenbasis))).T, fmt="% .6f", delimiter=" ")

  ###################################################################################
  # Initialize the density matrix, Rho

  Rho = np.zeros((dim, dim), dtype="complex128")
  Rho[0,0] = 1.0 

  ###################################################################################
  #

  D_level= np.zeros((nLevels, params.spinDegeneracy, dim, dim))
  for i in range(nLevels):
    for sigma in range(params.spinDegeneracy):
      D_level[i, sigma] = np.einsum('ba, bc, cd->ad', evecs, d_level[i, sigma], evecs)

  D_level_L = np.zeros((nLevels, params.spinDegeneracy, dim, dim))
  D_level_Ltilde = np.zeros((nLevels, params.spinDegeneracy, dim, dim))
  D_level_R = np.zeros((nLevels, params.spinDegeneracy, dim, dim))
  D_level_Rtilde = np.zeros((nLevels, params.spinDegeneracy, dim, dim))

  ###################################################################################
  # Main computational part of the program

  # Current and total population arrays to be filled in during propogation
  I_current = np.zeros(bias.shape[0])
  Spin_total = np.zeros(bias.shape[0])
  Pop_total = np.zeros(bias.shape[0])
  current_time_convergence = np.zeros(params.nTimeSteps)

  for i_bias in range(bias.shape[0]): 
    mL = params.biasEnergyCenter 
    mR = ( bias[i_bias] + params.biasEnergyCenter )
    for i_level in range(nLevels): 
      for sigma in range(params.spinDegeneracy): 
        for i in range(dim):
          for j in range(dim):
            D_level_Ltilde[i_level, sigma, i ,j] = D_level[i_level, sigma, i, j]*(1.0-fermi(evals[j]-evals[i], mL, params.kT))
            D_level_L[i_level, sigma, i, j] = D_level[i_level, sigma, i, j]*fermi(evals[j]-evals[i], mL, params.kT)
     
            D_level_Rtilde[i_level, sigma, i, j] = D_level[i_level, sigma, i, j]*(1.0-fermi(evals[j]-evals[i], mR, params.kT))
            D_level_R[i_level, sigma, i, j] = D_level[i_level, sigma, i, j]*fermi(evals[j]-evals[i], mR, params.kT)
      
    for i_level in range(nLevels): 
      for sigma in range(params.spinDegeneracy): 
        D_level_Ltilde[i_level, sigma] = np.einsum('ab,bc,dc->ad', evecs,D_level_Ltilde[i_level, sigma], evecs)
        D_level_L[i_level, sigma]= np.einsum('ab,bc,dc->ad', evecs, D_level_L[i_level, sigma], evecs)
        D_level_R[i_level, sigma]= np.einsum('ab,bc,dc->ad', evecs, D_level_R[i_level, sigma], evecs)
        D_level_Rtilde[i_level, sigma] = np.einsum('ab,bc,dc->ad', evecs, D_level_Rtilde[i_level, sigma], evecs)
            
    # Perform the time integration
    for t in range(params.nTimeSteps): 
      
      # Propogate the system density matrix, Rho, using the 4th-order Runge Kutta method
      Rho = lind_python2.RK4(D_level_Ltilde, D_level_L, D_level_Rtilde, D_level_R, d_level, Ham, Rho, params.dt, gammaLeft, gammaRight)

      if (i_bias == (bias.shape[0] / 2 )):
      	current_time_convergence[t] = calcCurrent(D_level_Ltilde, D_level_L, d_level, num_operator_total, Ham, Rho, gammaLeft)

    # Make sure the current is converged with respect to the time propogation length
    if (i_bias == (bias.shape[0] / 2 )):
      np.savetxt(str(bias[i_bias]) + "_currentTimeConvergence.dat", current_time_convergence, fmt="% .8f", delimiter=" ")

    # Calculate the steady-state current
    I_current[i_bias] = calcCurrent(D_level_Ltilde, D_level_L, d_level, num_operator_total, Ham, Rho, gammaLeft)

    # Store the total population
    Pop_total[i_bias] = np.real(np.einsum('ab,ba->', Rho, num_operator_total))

    # Calculate the total spin, S^2
    Spin_total[i_bias] = np.real(np.einsum('ab,ba->', Rho, S2_operator_total))

    # Print the steady-state solution of Rho
    np.savetxt(str(bias[i_bias]) + "_steadyStateDensityMatrix.dat", Rho, fmt="% .6f", delimiter=" ")
    np.savetxt(str(bias[i_bias]) + "_eigenSteadyStateDensityMatrix.dat", np.einsum('ba,bc,cd->ad', evecs, Rho, evecs), fmt="% .6f", delimiter=" ")

  ###################################################################################
  # Print the final results

  np.savetxt("currentVersusBias.dat", np.vstack((bias, I_current)).T, fmt="% .8f", delimiter=" ")
  np.savetxt("populationVersusBias.dat", np.vstack((bias, Pop_total)).T, fmt="% .8f", delimiter=" ")
  np.savetxt("spinVersusBias.dat", np.vstack((bias, Spin_total)).T, fmt="% .8f", delimiter=" ")

  ###################################################################################

###################################################################################

if __name__ == "__main__":
  main()

###################################################################################
