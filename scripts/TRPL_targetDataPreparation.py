###################################################################################
"""
	
	Program that performs the target data preparation for an optimization of the 
	kinetic rates and initial conditions of a classical master equation (i.e. a Markov
	chain Monte Carlo simulation). 

	For now, the file is specific to taking integrated time-resolved photoluminescence 
	data stored in (possibly multiple) comma separated csv files in which the 
	first column contains the times in picoseconds and the values in the other
	columns are the raw photon counts.  

	Author: John P. Philbin
	Last Modified: March 9th 2020

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
		--targetDataGlobalScaling=1.0 (-10.0, 30.0, ...)

	"""
	parser = ArgumentParser(fromfile_prefix_chars="@")

	parser.add_argument("--rawDataFilenameList", nargs="+", action="append", type=str,
							help="a list giving the amount to scale all target data values by")

	parser.add_argument("--nInitialConditions", type=int, default=1, 
							help="an into for the number of initial conditions (i.e. number of csv files)")

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
	entireRawData = np.genfromtxt(rawDataFilenameList[0][0], dtype=float, delimiter=",", skip_footer=0, 
								  missing_values=["NaN"], filling_values=[0.0])
	entireRawData = np.nan_to_num(entireRawData)

	# First column of the raw data will contain the times in picoseconds
	rawDataTimes = np.copy(entireRawData[:, 0])

	# All rows and columns excluding the first column 
	rawData = np.delete(np.copy(entireRawData), 0, axis=1)

	return np.array(entireRawData, ndmin=1), np.array(rawData, ndmin=1), np.array(rawDataTimes, ndmin=1)

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

	entireRawData, rawData, rawDataTimes = readRawDataFiles(params.rawDataFilenameList)
	print(entireRawData.shape)
	print(rawData.shape)
	print(rawDataTimes.shape)


	###############################################################################
	# Scale the photoluminescence intensity into <N> 

	targetData = ( params.targetDataGlobalScaling * rawData )

	np.savetxt("TARGETDATA.dat", targetData, fmt="% .6f", delimiter=" ")
	np.savetxt("TIMES.dat", rawDataTimes, fmt="% .6f", delimiter=" ")

	###############################################################################
	# Determine maximum and long time signals for each initial condition

	#maximumSignal = np.amax(np.abs(targetData), axis=1)
	#indexOfMaxSignal = np.argmax(np.abs(targetData), axis=1)
	#timeOfMaxSignal = np.zeros(indexOfMaxSignal.shape)
	#longTimeSignal = np.average(np.abs(targetData[:,-5:]), axis=1)
	#decayRatio = np.divide(maximumSignal, longTimeSignal)
	#for i in range(0, timeOfMaxSignal.shape[0]):
	#	timeOfMaxSignal[i] = rawDataTimes[indexOfMaxSignal[i]]
	#print("Time of maximum signal:\n")
	#print(timeOfMaxSignal)
	#print("Maximum signal:\n")	
	#print(maximumSignal)
	#print("Long time signal:\n")
	#print(longTimeSignal)
	#print("Long time decay ratio:\n")
	#print(decayRatio)

	###############################################################################
	# Write scaled alltargetData.par file

	#allTargetData = np.concatenate((np.reshape(np.transpose(rawDataTimes[0, :]), (rawDataTimes.shape[1], 1)), 
	#									np.transpose(targetData)), axis=1)
	#np.savetxt("allTargetData.par", allTargetData, fmt="% .8f", delimiter=" ")

	###############################################################################
	# Remove times below minimum and maximum times values and write final targetData.par file

	#finalTargetData = np.delete(allTargetData, np.where(allTargetData[:, 0] < params.minTime), axis=0)
	#finalTargetData = np.delete(finalTargetData, np.where(finalTargetData[:, 0] > params.maxTime), axis=0)
	#np.savetxt("targetData.par", finalTargetData, fmt="% .8f", delimiter=" ")

###################################################################################

if __name__ == "__main__":
	main()

###################################################################################
