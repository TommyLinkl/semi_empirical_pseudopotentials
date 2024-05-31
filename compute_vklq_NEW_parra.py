'''
Dipti Jasrasaria
September 2020

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# IMPORT MODULES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''

import os
import numpy as np
import argparse
import multiprocessing as mp
from functools import partial
'''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CONSTANTS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''
pstoau = 41341.373335182114 # 1ps in atomic time units

autoev = 27.211324570273 # 1 atomic unit (Hartree) in eV
autoj = 4.3597447222071e-18 # 1 Hartree in J
evtoj = 1.602176634e-19 # 1 eV in J

autom = 5.291772083e-11 # 1 atomic distance unit in meters

kb_ev = 8.617333262145e-5 # Boltzmann constant (eV*K)
kb_au = kb_ev/autoev # Boltzmann constant (Hartree*K)
kb_j = kb_ev*evtoj # Boltzmann constant (J*K)

hbar = 1.054571817e-34 # J*s 

'''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# COMMAND LINE ARGUMENTS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''

def buildParser():

    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--nAtom', type=int, default='0',
            help='number of semiconductor atoms in quantum dot structure')
    parser.add_argument('--nElec', type=int, default='0',
            help='number of electrons used in BSE')
    parser.add_argument('--nHole', type=int, default='0',
            help='number of holes used in BSE')
    parser.add_argument('--nExc', type=int, default='0',
            help='number of excitonic states to model')
    parser.add_argument('--temp', type=float, default='300',
            help='temperature (K)')
    return parser

'''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# READ INPUT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# read in BSE coefficients from intTwoParticleCoefficients.dat

def read_coeffs(nElec, nHole, nExc):
    filename = 'intTwoParticleCoefficients.dat'
    if not os.path.exists(filename):
        print('Cannot find ' + filename + '! Exiting...\n')
        exit()

    c = np.array([line.strip().split()[1] for line in open(filename, 'r')])
    if c.shape[0] != (nElec*nHole)**2:
        print('Error reading ' + filename + '! Exiting...\n')
        exit()
        
    c = (c[:nExc*nElec*nHole]).astype(np.float)
    return c

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # read masses of semiconductor atoms from conf.par

def read_mass(natoms):
    filename = 'conf.par'
    if not os.path.exists(filename):
        print('Cannot find ' + filename + '! Exiting...\n')
        exit()

    lines = np.array([line.strip().split() for line in open(filename, 'r')])[1:]
    at_type = np.array([line[0] for line in lines])[:natoms]
    
    # make sure no passivation atoms
    if ('P1' == at_type).any() or ('P2' == at_type).any():
        print('Error reading ' + filename + '! Exiting...\n')
        exit()

    mass = np.empty(natoms, dtype=np.float)
    cd_idx = np.where(at_type == 'Cd')[0]
    se_idx = np.where(at_type == 'Se')[0]
    s_idx = np.where(at_type == 'S')[0]
    if (cd_idx.shape[0] + se_idx.shape[0] + s_idx.shape[0] != natoms):
        print('Error reading ' + filename + '! Exiting...\n')
        exit()
    mass[cd_idx] = 112.411
    mass[se_idx] = 78.960
    mass[s_idx] = 32.065
    mass *= (1e-3/6.022140857e23) # kg

    return mass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # read in phonon modes from eig.dat
# i.e. eigenvectors of dynamical matrix to transform to phonon mode coordinates
# q = U*r

def read_phononModes(natoms):
    filename = 'eig.dat'
    if not os.path.exists(filename):
        print('Cannot find ' + filename + '! Exiting...\n')
        exit()

    pmodes = np.array([line.strip().split() for line in open(filename, 'r')]).astype(np.float) #.transpose()
    if (pmodes.shape[0] != 3*natoms) or (pmodes.shape[1] != 3*natoms):
        print ('Error reading ' + filename + '! Exiting...\n')
        exit()

    return pmodes

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # read in nonadiabatic coupling matrix elements from Vab-diabatic.dat or Vij-diabatic.dat

def read_DiabCoupling(natoms, nQPs, Vrs_file):

    if not os.path.exists(Vrs_file):
        print('Cannot find ' + Vrs_file + '! Exiting...\n')
        exit()

    lines = np.array([line.strip().split() for line in open(Vrs_file, 'r')])
    qp1Idx = np.array([line[6] for line in lines]).astype(np.int)
    qp2Idx = np.array([line[8] for line in lines]).astype(np.int)
    lines = lines[np.where((qp1Idx < nQPs) & (qp2Idx < nQPs))[0]][:nQPs*nQPs*natoms]
    if len(lines) != (natoms*nQPs*nQPs):
        print('Error reading in ' + Vrs_file + '! Exiting...\n')
        exit()
    
    Vrsx = np.array([line[10] for line in lines]).astype(np.float).reshape((natoms, nQPs, nQPs))
    Vrsy = np.array([line[11] for line in lines]).astype(np.float).reshape((natoms, nQPs, nQPs))
    Vrsz = np.array([line[12] for line in lines]).astype(np.float).reshape((natoms, nQPs, nQPs))
    Vrs = np.zeros((3*natoms, nQPs, nQPs))
    Vrs[::3,:,:] = Vrsx
    Vrs[1::3,:,:] = Vrsy
    Vrs[2::3,:,:] = Vrsz
    return Vrs

'''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# COMPUTE EXCITON-PHONON COUPLING
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # compute electron contribution of
# exciton-exciton diabatic coupling in atomic coordinates (Vkle)
def compute_Vkle(idx, natoms, nElec, nHole, coeffs, Vab):

    V = np.zeros(3*natoms)
    kk = idx[0]*nElec*nHole
    ll = idx[1]*nElec*nHole

    for i in range(nHole):
        cai = coeffs[kk+i*nElec:kk+(i+1)*nElec]
        cbi = coeffs[ll+i*nElec:ll+(i+1)*nElec]
        V += np.sum(np.einsum('i,j->ij', cai, cbi)*Vab, axis=(1,2))

    return (idx, V)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # compute hole contribution of
# exciton-exciton diabatic coupling in atomic coordinates (Vklh)
def compute_Vklh(idx, natoms, nElec, nHole, coeffs, Vij):

    V = np.zeros(3*natoms)
    kk = idx[0]*nElec*nHole
    kk1 = (idx[0]+1)*nElec*nHole
    ll = idx[1]*nElec*nHole
    ll1 = (idx[1]+1)*nElec*nHole

    for a in range(nElec):
        cai = coeffs[kk:kk1][a::nElec]
        caj = coeffs[ll:ll1][a::nElec]
        V += np.sum(np.einsum('i,j->ij', cai, caj)*Vij, axis=(1,2))

    return (idx, V)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # compute exciton-exciton diabatic coupling in atomic coordinates
# transform and sum single-particle couplings to exitonic basis

def compute_excCouplingAtomic(natoms, nElec, nHole, nExc):

    # read in Vkl from file
    if os.path.exists('Vkl-diabatic.dat'):
        Vkl = np.array([line.strip().split()[-1] for line in open('Vkl-diabatic.dat', 'r')]).astype(np.float)
        # check if excitonic diabatic couplings computed for nExc
        idx = np.sqrt(Vkl.shape[0]/(3*natoms)).astype(np.int)
        if idx >= nExc:
            Vkl = Vkl.reshape((3*natoms, idx, idx))
            return Vkl
        else:
            Vkle = np.zeros((3*natoms, nExc, nExc))
            tmp = np.array([line.strip().split()[-1] for line in open('Vkle-diabatic.dat', 'r')]).astype(np.float)
            Vkle[:,:idx,:idx] = tmp.reshape((3*natoms, idx, idx))

            Vklh = np.zeros((3*natoms, nExc, nExc))
            tmp = np.array([line.strip().split()[-1] for line in open('Vklh-diabatic.dat', 'r')]).astype(np.float)
            Vklh[:,:idx,:idx] = tmp.reshape((3*natoms, idx, idx))
    # compute Vkl
    else:
        idx = 0
        Vkle = np.zeros((3*natoms, nExc, nExc))
        Vklh = np.zeros((3*natoms, nExc, nExc))

    # BS coefficients
    c = read_coeffs(nElec, nHole, nExc)
    # quasiparticle diabatic couplings
    Vab = read_DiabCoupling(natoms, nElec, 'Vab-diabatic.dat')
    Vij = read_DiabCoupling(natoms, nHole, 'Vij-diabatic.dat')
    
#####here is paraller part begin Vkle####
    #get_indices
    idxs = [(k,l) for k in range(nExc) for l in range(k+1,nExc)] # upper triangular
    idxs_d = [(k,k) for k in range(nExc)] # diagonal

    pool = mp.Pool(mp.cpu_count())
    mats = pool.map(partial(compute_Vkle, natoms=natoms, nElec=nElec, nHole=nHole, coeffs=c, Vab=Vab), idxs)
    mats_d = pool.map(partial(compute_Vkle, natoms=natoms, nElec=nElec, nHole=nHole, coeffs=c, Vab=Vab), idxs_d)

    Vkle = np.empty((3*natoms, nExc, nExc))
    for i in range(len(mats)):
        (k,l) = mats[i][0]
        Vkle[:,k,l] = mats[i][1]
        Vkle[:,l,k] = mats[i][1]
    for i in range(len(mats_d)):
        k = mats_d[i][0][0]
        Vkle[:,k,k] = mats_d[i][1]

    mats = pool.map(partial(compute_Vklh, natoms=natoms, nElec=nElec, nHole=nHole, coeffs=c, Vij=Vij), idxs)
    mats_d = pool.map(partial(compute_Vklh, natoms=natoms, nElec=nElec, nHole=nHole, coeffs=c, Vij=Vij), idxs_d)

    Vklh = np.empty((3*natoms, nExc, nExc))
    for i in range(len(mats)):
        (k,l) = mats[i][0]
        Vklh[:,k,l] = mats[i][1]
        Vklh[:,l,k] = mats[i][1]
    for i in range(len(mats_d)):
        k = mats_d[i][0][0]
        Vklh[:,k,k] = mats_d[i][1]

    pool.close()
    pool.join()

    # electron contribution
    #for k in range(idx, nExc):
    #   kk = k*nElec*nHole
    #   for l in range(idx, nExc):
    #       ll = l*nElec*nHole
    #       for i in range(nHole):
    #           cai = c[kk+i*nElec:kk+(i+1)*nElec]
    #           cbi = c[ll+i*nElec:ll+(i+1)*nElec]
    #           Vkle[:,k,l] += np.sum(np.einsum('i,j->ij', cai, cbi)*Vab, axis=(1,2))
    #with open('Vkle-diabatic.dat', 'w') as f:
    #    for n in range(Vkle.shape[0]):
    #        for i in range(Vkle.shape[1]):
    #            for j in range(Vkle.shape[2]):
    #                f.write(str(n) + ' ' + str(i) + ' ' + str(j) + ' ' + str(Vkle[n,i,j]) + '\n')

    # hole contribution
    # for k in range(idx,nExc):
    #   kk = k*nElec*nHole
    #   kk1 = (k+1)*nElec*nHole
    #   for l in range(idx,nExc):
    #       ll = l*nElec*nHole
    #       ll1 = (l+1)*nElec*nHole
    #       for a in range(nElec):
    #           cai = c[kk:kk1][a::nElec]
    #           caj = c[ll:ll1][a::nElec]
    #           Vklh[:,k,l] += np.sum(np.einsum('i,j->ij', cai, caj)*Vij, axis=(1,2))
    #with open('Vklh-diabatic.dat', 'w') as f:
    #    for n in range(Vklh.shape[0]):
    #        for i in range(Vklh.shape[1]):
    #            for j in range(Vklh.shape[2]):
    #                f.write(str(n) + ' ' + str(i) + ' ' + str(j) + ' ' + str(Vklh[n,i,j]) + '\n')

    Vkl = Vkle-Vklh 
    Vkl *= autoj/autom # J / m
    with open('Vkl-diabatic.dat', 'w') as f:
        for n in range(Vkl.shape[0]):
            for i in range(Vkl.shape[1]):
                for j in range(Vkl.shape[2]):
                    f.write(str(n) + ' ' + str(i) + ' ' + str(j) + ' ' + str(Vkl[n,i,j]) + '\n')

    return Vkl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # compute exciton-exciton diabatic coupling in phonon coordinates
# transform from atomic coordinates to phonon mode coordinates

def compute_excCouplingPhonon(natoms, nElec, nHole, nExc):

    # read in Vkl from file
    if os.path.exists('Vklq-diabatic.dat'):
        Vklq = np.array([line.strip().split()[-1] for line in open('Vklq-diabatic.dat', 'r')]).astype(np.float)
        # check if excitonic diabatic couplings computed for nExc
        idx = np.sqrt(Vklq.shape[0]/(3*natoms)).astype(np.int)
        if idx >= nExc:
            Vklq = Vklq.reshape((3*natoms, idx, idx))
            return Vklq

    # Vkl in atomic coordinates ( J / m)
    Vkl = compute_excCouplingAtomic(natoms, nElec, nHole, nExc)

    # phonon modes
    pmodes = read_phononModes(natoms)
    pmodes_inv = pmodes.transpose()
    sr_mass3 = np.sqrt(np.repeat(read_mass(natoms), 3)) # sqrt(kg)

    # compute Vklq 
    Vklq = np.empty((3*natoms, nExc, nExc))
    for k in range(Vkl.shape[1]):
        for l in range(Vkl.shape[2]):
            Vklq[:,k,l] = np.matmul(pmodes_inv, Vkl[:,k,l]/sr_mass3) # J / sqrt(kg)*m

    with open('Vklq-diabatic.dat', 'w') as f:
        for n in range(Vklq.shape[0]):
            for i in range(Vklq.shape[1]):
                for j in range(Vklq.shape[2]):
                    f.write(str(n) + ' ' + str(i) + ' ' + str(j) + ' ' + str(Vklq[n,i,j]) + '\n')

    return Vklq

'''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MAIN
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''

# # # # parameters
parser = buildParser()
params = parser.parse_args()

vklq = compute_excCouplingPhonon(params.nAtom, params.nElec, params.nHole, params.nExc)[6:,:,:] # J / sqrt(kg)*m
