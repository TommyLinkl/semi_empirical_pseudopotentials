#######################################################################################
"""
	
	Script that creates kernel_mtxels.h5 file that can be used as input for the
	calculation of fully-correlated excitons, trions and biexcitons using 
	Berkeley GW in collaboration with Dr. Felipe H. da Jornada of the Louie Group

	Author: John P. Philbin
	Last Modified: August 28th, 2019

"""
#######################################################################################

# External libraries
import os
import time
import math
import h5py
import numpy as np 
from argparse import ArgumentParser

#######################################################################################

# Unit conversion
AUTOEV = np.float_(27.2114)

#######################################################################################

def buildParser():
	"""
	Parser for command line arguments
	Options:
		--mode: 
	"""
	parser = ArgumentParser()

	parser.add_argument("--calcType", type=str, default="negativeTrion", 
						help="allowed calcTypes are negativeTrion, positiveTrion and biexciton")
	parser.add_argument("--inputDirNum", type=str, default="1",	
						help="inputDirNum (str) to read input files from (ex: 2 for ./2_elecElec/")
	parser.add_argument("--nElecs", type=int, default=1, help="number of electron states")
	parser.add_argument("--nHoles", type=int, default=1, help="number of hole states")
	parser.add_argument("--nKPoints", type=int, default=1, help="number of k-points")
	parser.add_argument("--outputDirNum", type=str, default="0", 
						help="outputDirNum (str) to put output file in (ex: 2 for ./2_trions/")

	return parser

#######################################################################################

def readQuasiparticleEnergies(params):
	"""
	
	"""

	qpEnergies = np.zeros((params.nKPoints, params.nElecs+params.nHoles), dtype=float)

	inputFileName = params.inputDirNum + "_elecHole/holeElecEigenstates.dat"
	with open(inputFileName) as f:
		iQP = 0
		for line in f:
			row = line.split()
			if (iQP < params.nHoles):
				qpEnergies[0, params.nHoles-iQP-1] = np.float_(row[2])*AUTOEV
			else:
				qpEnergies[0, iQP] = np.float_(row[2])*AUTOEV
			iQP += 1

	return qpEnergies

#######################################################################################

def readCoulombMtxels(params, calcType):
	"""

	"""
	
	# Determine size of mtxels and inputFileName based on calcType and input params
	nk = params.nKPoints
	if calcType == "W_cv" or calcType == "v_cv":
		if calcType == "W_cv":
			inputFileName = params.inputDirNum + "_elecHole/subset_Wrsut.dat"
		else:
			inputFileName = params.inputDirNum + "_elecHole/subset_Vrsut.dat"
		nb1 = nb3 = params.nElecs
		nb2 = nb4 = params.nHoles
		ib1 = 1 # ic1,  elec 1
		ib2 = 0 # iv1,  hole 1
		ib3 = 3 # ic1', elec 1'
		ib4 = 2 # iv1', hole 1'
		iValue = 8
	elif calcType == "W_cc": # <c1c2|W|c1'c2'>
		inputFileName = params.inputDirNum + "_elecElec/Wrsut.dat"
		nb1 = nb2 = nb3 = nb4 = params.nElecs
		ib1 = 0 # ic1,  elec 1 
		ib2 = 1 # ic2,  elec 2
		ib3 = 2 # ic1', elec 1'
		ib4 = 3 # ic2', elec 2'
		iValue = 4
	elif calcType == "W_vv": 
		inputFileName = params.inputDirNum + "_holeHole/Wrsut.dat"
		nb1 = nb2 = nb3 = nb4 = params.nHoles
		ib1 = 0 # ic1,  hole 1 
		ib2 = 1 # ic2,  hole 2
		ib3 = 2 # ic1', hole 1'
		ib4 = 3 # ic2', hole 2'
		iValue = 4
	else:
		print "Error in calling readCoulombMtxels, invalid calcType"
		return 0

	mtxels = np.zeros((2, nk, nb1, nb2, nk, nb3, nb4, nk), dtype=float)
	# mtxels[1,:] = 0.0 as this is the imaginary part

	with open(inputFileName) as f:
		for line in f:
			row = line.split()
			mtxels[0, 0, int(row[ib1]), int(row[ib2]), 0, int(row[ib3]), int(row[ib4]), 0] = 2.0*np.float_(row[iValue])

	# Transpose to get the correct ordering for Felipe's code
	mtxels = np.transpose(mtxels)

	return mtxels

#######################################################################################

def printH5PYDatasetInfo(name, node):
	"""
	
	"""

	if isinstance(node, h5py.Dataset):
		print(str(node.name) + " with shape " + str(node.shape))

	return

#######################################################################################

def main():
	"""
	
	"""

	###################################################################################
	# Parse command line arguments
	parser = buildParser()
	params = parser.parse_args()

	###################################################################################
	# Create the hdf5 file for which the results will be stored

	if params.outputDirNum == "0":
		outputDir = "./" + params.inputDirNum + "_trions/"
	else:
		outputDir = "./" + params.outputDirNum + "_trions/"
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	finalFile = h5py.File(outputDir + "kernel_mtxels.h5", "w")

	###################################################################################
	# Create integer datasets based off command line input arguments

	# dataset for the number of conduction bands / electron states 
	nElecs = np.int_(params.nElecs).reshape(1)
	nc = finalFile.create_dataset("nc", data=nElecs)

	# dataset for the number of valence bands / hole states 
	nHoles = np.int_(params.nHoles).reshape(1)
	nv = finalFile.create_dataset("nv", data=nHoles)

	# dataset for the number of k-points (should always be 1)
	nKPoints = np.int_(params.nKPoints).reshape(1)
	nk = finalFile.create_dataset("nk", data=nKPoints)

	###################################################################################
	# Create double (float_) dataset for the k-points under Gamma point only assumption

	# dataset for the k-points 
	kPoints = np.float_([0.0, 0.0, 0.0]).reshape(1,3)
	kpts = finalFile.create_dataset("kpts", data=kPoints) 

	###################################################################################
	# Create double (float_) datasets that are read in from files

	# array and dataset for the quasiparticle energies 
	qpEnergies = readQuasiparticleEnergies(params)
	eqp = finalFile.create_dataset("eqp", data=qpEnergies) 

	# array and dataset for the screened electron-electron Coulomb matrix elements
	elecElecMtxels = readCoulombMtxels(params, "W_cc")
	W_cc = finalFile.create_dataset("W_cc", data=elecElecMtxels)

	# array and dataset for the screened electron-hole Coulomb matrix elements
	elecHoleMtxels = readCoulombMtxels(params, "W_cv")
	W_cv = finalFile.create_dataset("W_cv", data=elecHoleMtxels)

	# array and dataset for the bare electron-hole Coulomb matrix elements
	bareElecHoleMtxels = readCoulombMtxels(params, "v_cv")
	v_cv = finalFile.create_dataset("v_cv", data=bareElecHoleMtxels)

	# array and dataset for the screened hole-hole Coulomb matrix elements
	holeHoleMtxels = readCoulombMtxels(params, "W_vv")
	W_vv = finalFile.create_dataset("W_vv", data=holeHoleMtxels)

	###################################################################################

	# Print summary of input and output file
	print "\nCreated the following file:"
	print outputDir + "kernel_mtxels.h5\n"
	print "Containing the following datasets:"
	finalFile.visititems(printH5PYDatasetInfo)
	print ""

	###################################################################################

if __name__ == "__main__":
	main()

#######################################################################################
