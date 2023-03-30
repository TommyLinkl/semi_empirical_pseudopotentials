###################################################################################
"""
	
	Program that performs the target data preparation for an optimization of the 
	kinetic rates and initial conditions of a classical master equation (i.e. a Markov
	chain Monte Carlo simulation). 

	For now, the file is specific to taking transient absorption data stored in 
	(possibly multiple) comma separated csv files in which the first column contains 
	the detection wavelengths in nm, the first row contains the times in picoseconds
	and the values in the matrix are the raw transient absorption signals.  

	Author: John P. Philbin
	Last Modified: December 11th 2019

	Input files:
		dataPrepInput.par -> one argument per line, example: --outputDirNum=1

	Output files:
		targetData.par ->

	Example run commands:
		$ conda activate
		$ python targetDataPreparation.py @dataPrepInput.par

"""
###################################################################################

# External libraries
import os
import math
import shutil
import numpy as np 
from argparse import ArgumentParser

###################################################################################

# Functions

def buildParser():
	"""
	Parser for command line arguments or from an input file (e.g. dataPrepInput.par)
	Options:
		--nInitialConditions=1 (2, 3, ...)
		--nTargetPopulations=1 (2, 3, ...)
		--targetDataGlobalScaling=1.0 (-10.0, 30.0, ...)

	"""
	parser = ArgumentParser(fromfile_prefix_chars="@")

	parser.add_argument("--rawDataFilenameList", nargs="+", action="append", type=str,
							help="a list giving the amount to scale all target data values by")

	parser.add_argument("--nInitialConditions", type=int, default=1, 
							help="an into for the number of initial conditions (i.e. number of csv files)")
	parser.add_argument("--nTargetPopulations", type=int, default=1, 
							help="an int for the number of the target populations that we want to fit to")

	parser.add_argument("--minWavelength", type=float, default=499.0, 
							help="minimum wavelength for which averaging will be done over")
	parser.add_argument("--maxWavelength", type=float, default=501.0, 
							help="maximum wavelength for which averaging will be done over")	
	parser.add_argument("--minTime", type=float, default=(-1.0*float("inf")), 
							help="minimum time that will be included in the final targetData.par")
	parser.add_argument("--maxTime", type=float, default=float("inf"), 
							help="maximum time that will be included in the final targetData.par")							
	parser.add_argument("--targetDataGlobalScaling", type=float, default=1.0, 
							help="a float giving the amount to scale all target data values by")	

	return parser

###################################################################################

def readRawDataFiles(rawDataFilenameList):
	""" 
	Input:
		- rawDataFilenameList = the name of the files that contains the raw
	Returns:
		- entireRawData = the complete csv files
		- targetData = a numpy array of shape (nTimeSteps, nInitialConditions) that
					   contains the target time evolution of the target population
		- targetDataTimes = a numpy array of shate (nTimeSteps) that contains the times
							at which the data was collected
		- targetDataTimes = a numpy array of shate (nTimeSteps) that contains the times
							at which the data was collected
	"""
	
	# Read in the entire raw data csv file 
	for i, fileName in enumerate(rawDataFilenameList):
		tmpArray =  np.genfromtxt(fileName[0], dtype=float, delimiter=",", skip_footer=15, 
								  missing_values=["NaN"], filling_values=[0.0])
		if i == 0:
			entireRawData = np.empty((len(rawDataFilenameList), tmpArray.shape[0], tmpArray.shape[1]), dtype=float)
		entireRawData[i, :, :] = np.copy(tmpArray)
	entireRawData = np.nan_to_num(entireRawData)

	# First row of the raw data will contain the times in picoseconds
	rawDataTimes = np.copy(entireRawData[:, 0, 1:])

	# First column of the raw data will contain the detection wavelengths in nm 
	rawDataWavelengths = np.copy(entireRawData[:, 1:, 0])

	# All rows and columns excluding the first row and column 
	rawData = np.delete(np.copy(entireRawData), 0, axis=2)
	rawData = np.delete(rawData, 0, axis=1)

	return np.array(entireRawData, ndmin=1), np.array(rawData, ndmin=1), np.array(rawDataTimes, ndmin=1), np.array(rawDataWavelengths, ndmin=1)

###################################################################################

def main():	
	""" 

	"""
	
	###############################################################################
	# Parse command line arguments

	parser = buildParser()
	params = parser.parse_args()

	###############################################################################
	# Reads in the raw csv file(s)

	entireRawData, rawData, rawDataTimes, rawDataWavelengths = readRawDataFiles(params.rawDataFilenameList)

	###############################################################################
	# Extract and average over desired wavelength window for each initial condition

	targetData = np.zeros_like(rawDataTimes)

	for i in range(0, rawData.shape[0]):
		nWavelengthsInRange = 0.0
		for j in range(0, rawData.shape[1]):
			if (rawDataWavelengths[i, j] > params.minWavelength and rawDataWavelengths[i, j] < params.maxWavelength):
				targetData[i, :] = np.add(targetData[i, :], rawData[i, j, :])
				nWavelengthsInRange += 1.0
		if (nWavelengthsInRange > 0.5):
			targetData[i, :] /= nWavelengthsInRange
	targetData *= params.targetDataGlobalScaling

	###############################################################################
	# Determine maximum and long time signals for each initial condition

	maximumSignal = np.amax(np.abs(targetData), axis=1)
	indexOfMaxSignal = np.argmax(np.abs(targetData), axis=1)
	timeOfMaxSignal = np.zeros(indexOfMaxSignal.shape)
	longTimeSignal = np.average(np.abs(targetData[:,-5:]), axis=1)
	decayRatio = np.divide(maximumSignal, longTimeSignal)
	for i in range(0, timeOfMaxSignal.shape[0]):
		timeOfMaxSignal[i] = rawDataTimes[i, indexOfMaxSignal[i]]
	print(timeOfMaxSignal)
	print(maximumSignal)
	print(longTimeSignal)
	print(decayRatio)

	###############################################################################
	# Write scaled alltargetData.par file

	allTargetData = np.concatenate((np.reshape(np.transpose(rawDataTimes[0, :]), (rawDataTimes.shape[1], 1)), 
										np.transpose(targetData)), axis=1)
	np.savetxt("allTargetData.par", allTargetData, fmt="% .8f", delimiter=" ")

	###############################################################################
	# Remove times below minimum and maximum times values and write final targetData.par file

	finalTargetData = np.delete(allTargetData, np.where(allTargetData[:, 0] < params.minTime), axis=0)
	finalTargetData = np.delete(finalTargetData, np.where(finalTargetData[:, 0] > params.maxTime), axis=0)
	np.savetxt("targetData.par", finalTargetData, fmt="% .8f", delimiter=" ")

###################################################################################

if __name__ == "__main__":
	main()

###################################################################################
