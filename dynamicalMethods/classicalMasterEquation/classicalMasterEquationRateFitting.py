from __future__ import print_function

###################################################################################
"""
	
	Program that performs a Monte Carlo fit for the 
	rate constnats of a kinetic equation to target data 
	using a Monte Carlo scheme

	Author: John P. Philbin
	Last Modified: December 11th 2019

	Input files:
		input.par -> one argument per line, example: --outputDirNum=1
		params.par ->
		system.par ->
		initialConditions.par -> 
		rates.par ->
		targetData.par -> 

	Output files:
		xx ->

	Example run commands:
		$ conda activate
		$ python classicalMasterEquationRateFitting.py @input.par
		$ python classicalMasterEquationRateFitting.py --nRates=6 --nInitialCondition=2

"""
###################################################################################

# External libraries
import os
import sys
import subprocess
import time
import math
import shutil
import numpy as np 
from argparse import ArgumentParser
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import normalize
#from skopt import gp_minimize
#from skopt import plot_convergence

###################################################################################

# Functions

def buildParser():
	"""
	Parser for command line arguments
	Options:
		--mode: monteCarloFit
	"""
	parser = ArgumentParser(fromfile_prefix_chars="@")

	parser.add_argument("--searchMethod", type=str, default="monteCarlo",
							help="search method, either grid or monteCarlo")

	parser.add_argument("--outputDirNum", type=int, default=1, 
							help="outputDirNum (int) to put output file in (ex: 2 for ./2_fit/")
	parser.add_argument("--mode", type=str, default="monteCarloFit", 
							help="allowed modes are monteCarloFit and singleIteration")
	parser.add_argument("--nTimeSteps", type=int, default=100, 
							help="number of time steps for each initial condition")
	parser.add_argument("--nRates", type=int, default=1, 
							help="number of rate constants to fit")
	parser.add_argument("--nInitialConditions", type=int, default=1, 
							help="number of initial conditions (i.e. populations)")
	parser.add_argument("--iTargetPopulation", type=int, default=1, 
							help="index of the target population that we want to fit to")
	parser.add_argument("--nMonteCarloIterations", type=int, default=1, 
							help="number of Monte Carlo iterations (e.g. 0, 1, ... 1000, ...)")
	parser.add_argument("--monteCarloBeta", type=float, default=100.0, 
							help="beta used in deciding whether to keep or discard the new rate constants")
	parser.add_argument("--monteCarloStepSize", type=float, default=0.01, 
							help="size of the changes to the rate constants at each Monte Carlo iteration")
	parser.add_argument("--writeIterationStepFrequency", type=int, default=10, 
							help="append to iterations.dat every this number of steps")
	parser.add_argument("--targetDataGlobalScaling", type=float, default=1.0, 
							help="amount to scale all target data values by for the fitting")	
	parser.add_argument("--maxTargetDataGlobalScaling", type=float, default=1.0, 
							help="the maximum amount to scale all target data values by for the fitting")
	parser.add_argument("--minTargetDataGlobalScaling", type=float, default=1.0, 
							help="the minimum amount to scale all target data values by for the fitting")
	parser.add_argument("--targetDataGlobalScalingStepSize", type=float, default=0.0, 
							help="step size for grid search")
	parser.add_argument("--machineEpsilon", type=float, default=1.0e-10, 
							help="machine epsilon, number smallers than these are considered equal to 0.0")

	return parser

###################################################################################

def readFittingInputParams(fittingParamFileName):
	""" 
	Input:
		- fittingParamFileName = the name of the file that contains the initial guesses
							     for the base rates and by how much to change them each iteration 
	Returns:
		- fittingRateParams = a numpy array of shape (nRates) that
					   		  contains the initial guess for the initial populations
		- fittingStepSizeParams = a numpy array of shape (nRates) that
					   		  contains the step size for each base rate
	"""

	fittingRateParamsAll = np.loadtxt(fittingParamFileName, dtype=float)
	#fittingRateParamsAll = np.array(fittingRateParamsAll, ndmin=2)

	if fittingRateParamsAll.shape[1] == 2:
		fittingRateParams = np.loadtxt(fittingParamFileName, dtype=float, usecols=(0,))
		fittingStepSizeParams = np.loadtxt(fittingParamFileName, dtype=float, usecols=(1,))
	else: # fittingRateParamsAll.shape[1] == 3
		fittingRateParams = np.loadtxt(fittingParamFileName, dtype=float, usecols=(0,1))
		fittingStepSizeParams = np.loadtxt(fittingParamFileName, dtype=float, usecols=(2,))		

	return np.array(fittingRateParams, ndmin=1), np.array(fittingStepSizeParams, ndmin=1) 

###################################################################################

def readInitialConditionParams(initialConditionsParamFileName):
	""" 
	Input:
		- initialConditionsParamFileName = the name of the file that contains the initial guesses
							   for the initial conditions and by how much to change them each iteration 
	Returns:
		- initialConditions = a numpy array of shape (nInitialConditions) that contains
					  		  the initial guess for the initial populations
		- initialConditionsStepSize = a numpy array of shape (nInitialConditions) that contains
					  		  the step size for each initial populations
	"""

	initialConditionsAll = np.loadtxt(initialConditionsParamFileName, dtype=float)
	#initialConditionsAll = np.array(initialConditionsAll, ndmin=2)

	if initialConditionsAll.shape[1] == 2:
		initialConditions = np.loadtxt(initialConditionsParamFileName, dtype=float, usecols=(0,))
		initialConditionsStepSize = np.loadtxt(initialConditionsParamFileName, dtype=float, usecols=(1,))
	else: # initialConditionsAll.shape[1] == 3 for grid search
		initialConditions = np.loadtxt(initialConditionsParamFileName, dtype=float, usecols=(0,1))
		initialConditionsStepSize = np.loadtxt(initialConditionsParamFileName, dtype=float, usecols=(2,))

	return np.array(initialConditions, ndmin=1), np.array(initialConditionsStepSize, ndmin=1)

###################################################################################

def readTargetDataFile(targetDataFileName):
	""" 
	Input:
		- targetDataFileName = the name of the file that contains the target time evolution
							   of just the target population for all of the initial conditions.
							   The first column should have the times at which data exists
							   and the remaining columns should have the populations for each
							   initial condition with spaces separating each column 
	Returns:
		- targetData = a numpy array of shape (nTimeSteps, nInitialConditions) that
					   contains the target time evolution of the target population
		- targetDataTimes = a numpy array of shate (nTimeSteps) that contains the times
							at which the data was collected
	"""
	
	# Read in the entire targetDataFile 
	targetDataAndTimes = np.loadtxt(targetDataFileName, dtype=float)
	
	# First column of targetData will contain the times, t
	targetDataTimes = targetDataAndTimes[:, 0]

	# All columns after the first contain the values, f(t)
	targetData = np.delete(targetDataAndTimes, 0, axis=1)

	return np.array(targetData, ndmin=1), np.array(targetDataTimes, ndmin=1)

###################################################################################

def updateClassicalMasterEquationSystemParFile(rateParams, systemParFileName):
	""" 

	"""

	with open(systemParFileName, "r") as f:
		systemFile = f.readlines()

	i = 0
	indexList = []
	for lineIndex, line in enumerate(systemFile):
		line = line.strip()
		if line.strip():
			row = line.split()
			if ( row[0] == "baseRate" ):
				indexList.append(lineIndex) 
				i += 1

	for j in range(0, i):
		systemFile[indexList[j]] = "baseRate = %.8f\n" % rateParams[j]

	with open("system.dat", "w") as f:
		f.writelines( systemFile )

	return

###################################################################################

def updateClassicalMasterEquationParamsParFile(initialCondition, paramsFileName):
	""" 

	"""

	with open(paramsFileName, "r") as f:
		paramsFile = f.readlines()

	index = 0
	seedIndex = 0
	for lineIndex, line in enumerate(paramsFile):
		line = line.strip()
		if line.strip():
			row = line.split()
			if ( row[0] == "aveInitialPop"):
				index = lineIndex 
			if ( row[0] == "seed"):
				seedIndex = lineIndex

	paramsFile[index] = "aveInitialPop = %.4f\n" % initialCondition
	paramsFile[seedIndex] = "seed = %ld\n" % np.random.randint(1000, high=10000000)

	with open("params.dat", "w") as f:
		f.writelines( paramsFile )

	return

###################################################################################

def monteCarloUpdate(originalParams, stepSize, nParams):
	""" 

	"""

	randomArray = np.subtract(np.multiply(2.0, np.random.rand(nParams)), 1.0) 

	dParams = np.multiply(stepSize, randomArray)

	newParams = np.add(originalParams, dParams)

	return newParams

###################################################################################

def solveClassicalMasterEquation(pathToCalculationDirectory):
	""" 

	"""

	currentDirectory = os.getcwd()

	os.chdir(pathToCalculationDirectory)

	executableCommand = "/home/jphilbin/programs/semi_empirical_pseudopotentials/dynamicalMethods"
	executableCommand = executableCommand + "/classicalMasterEquation/kinMCvDt.x"

	subprocess.call(executableCommand)

	populations = np.loadtxt("trajectories.dat", dtype=float)

	os.chdir(currentDirectory)

	return populations


def solveClassicalMasterEquationMultipleInitialConditions(data, baseRate, initialConditions, 
														  initialConditionDirs, systemParFileName, paramsFileName, 
														  targetDataTimes, targetData, params):
	"""

	"""

	updateClassicalMasterEquationSystemParFile(baseRate, systemParFileName)

	for i in range (0, params.nInitialConditions):
		updateClassicalMasterEquationParamsParFile(initialConditions[i], paramsFileName)
		
		# Copy input files needed for solving the classical master equation
		shutil.copyfile("system.dat", initConditionDirs[i] + "system.par")
		shutil.copyfile("params.dat", initConditionDirs[i] + "params.par")
	
		# Calls c executable kinMCvDt.x 
		allPopulations = solveClassicalMasterEquation(initConditionDirs[i])
	
		# Interpolate data to target data time points and store target population in data
		data[:, i] = np.interp(targetDataTimes, allPopulations[:, 0], allPopulations[:, params.iTargetPopulation])
	
	# Calculate the mean-squared error using scikit-learn library function
	meanSquaredError = mean_squared_error(data, targetData)

	return meanSquaredError


###################################################################################

def writeIteration(i, mse, baseRate, initialCondition, iterationFileName):
	""" 

	"""

	with open(iterationFileName, "a+") as f:
		f.write("%5d  %.8f " % (i, mse))
		for i in range(0, baseRate.shape[0]):
			f.write("%.8f " % baseRate[i])
		for i in range(0, initialCondition.shape[0]):
			f.write(" %.4f  " % initialCondition[i])
		f.write("\n")

	return

###################################################################################

def writeBestRateParams(bestMSE, bestRates, bestInitialConditions, outputFileName):
	""" 

	"""
	with open(outputFileName, "w") as f:
		f.write("Best mean-squared error = %.12f\n" % bestMSE)
		for i in range(0, bestRates.shape[0]):
			f.write("baseRate %2d = %.8f\n" % (i, bestRates[i]))
		for i in range(0, bestInitialConditions.shape[0]):
			f.write("initCond %2d = %.8f\n" % (i, bestInitialConditions[i]))

	return

###################################################################################
 

def printCombinationsToFile(arr, fileName): 
    """
	
	function to make combinations that contain one element from each of the given arrays
    
    """

    # file to print the result
    file = open(fileName, 'w')

    # number of arrays 
    n = len(arr) 
  
    # to keep track of next element  
    # in each of the n arrays 
    indices = [0 for i in range(n)] 
  
    while (1): 
  
        # prcurrent combination 
        for i in range(n): 
            print(arr[i][indices[i]], end=" ", file=file) 
        print(file=file) 
  
        # find the rightmost array that has more 
        # elements left after the current element 
        # in that array 
        next = n - 1
        while (next >= 0 and 
              (indices[next] + 1 >= len(arr[next]))): 
            next-=1
  
        # no such array is found so no more 
        # combinations left 
        if (next < 0): 
            return
  
        # if found move to next element in that 
        # array 
        indices[next] += 1
  
        # for all arrays to the right of this 
        # array current index again points to 
        # first element 
        for i in range(next + 1, n): 
            indices[i] = 0

###################################################################################

def main():	
	""" 

	"""
	
	###############################################################################
	# Parse command line arguments

	parser = buildParser()
	params = parser.parse_args()

	###############################################################################
	# Reads in the input file parameters and target data file
	
	baseRate, stepSize = readFittingInputParams("rates.par")

	initialConditions, initialConditionStepSize = readInitialConditionParams("initialConditions.par")

	targetData, targetDataTimes = readTargetDataFile("targetData.par")
	normedTargetData = normalize(targetData, axis=0, norm="l2")

	systemParFileName = "system.par"

	paramsFileName = "params.par"

	###############################################################################
	# Create the folder and files for which the results will be stored
	
	outputDir = "./" + str(params.outputDirNum) + "_fit/"
	while(os.path.exists(outputDir)):
		params.outputDirNum += 1
		outputDir = "./" + str(params.outputDirNum) + "_fit/"

	if not os.path.exists(outputDir):
		os.makedirs(outputDir)

	if os.path.isfile("input.par"):
		shutil.copyfile("input.par", outputDir + "input.dat")

	outputFileName = outputDir + "fit.dat"
	iterationFileName = outputDir + "iterations.dat"
	with open(iterationFileName, "w") as f:
		f.write("#iIter     MSE   ")
		for i in range(0, params.nRates):
			f.write("    Rate   ")
		for i in range(0, params.nInitialConditions):
			f.write("  <N(0)> ")
		f.write("\n")

	###############################################################################
	# Creates the folders in which the trajectories for each initial condition 

	initConditionDirs = []
	for i in range(0, params.nInitialConditions):
		initConditionDirs.append(outputDir + "initCondition_" + str(i) + "/")
		if not os.path.exists(initConditionDirs[i]):
			os.makedirs(initConditionDirs[i])
		shutil.copy(paramsFileName, initConditionDirs[i])
		
	###############################################################################
	# Perform Monte Carlo fitting

	bestMSE, currentMSE = float("inf"), float("inf")
	bestRates, currentRates = np.copy(baseRate), np.copy(baseRate)
	bestInitialConditions, currentInitialConditions = np.copy(initialConditions), np.copy(initialConditions)

	data = np.zeros_like(targetData)

	if params.searchMethod == "monteCarlo":
		for i in range (0, params.nMonteCarloIterations):

			currentRates, currentInitialConditions = np.copy(baseRate), np.copy(initialConditions)
			
			if i:
				baseRate = monteCarloUpdate(baseRate, stepSize, params.nRates)
				initialConditions = monteCarloUpdate(initialConditions, initialConditionStepSize, params.nInitialConditions)
			
			updateClassicalMasterEquationSystemParFile(baseRate, systemParFileName)
			
			for j in range (0, params.nInitialConditions):

				updateClassicalMasterEquationParamsParFile(initialConditions[j], paramsFileName)

				# Copy input files needed for solving the classical master equation
				shutil.copyfile("system.dat", initConditionDirs[j] + "system.par")
				shutil.copyfile("params.dat", initConditionDirs[j] + "params.par")

				# Calls c executable kinMCvDt.x 
				allPopulations = solveClassicalMasterEquation(initConditionDirs[j])

				# Interpolate data to target data time points and store target population in data
				data[:, j] = np.interp(targetDataTimes, allPopulations[:, 0], allPopulations[:, params.iTargetPopulation])
			
			# Calculate the mean-squared error using scikit-learn library function
			#newMSE = mean_squared_error(normalize(data, axis=0, norm="l2"), normedTargetData)/params.nInitialConditions
			newMSE = mean_squared_error(data, targetData)

			# Write current mean-squared error
			if (i % params.writeIterationStepFrequency == 0):
				writeIteration(i, newMSE, baseRate, initialConditions, iterationFileName)

			# Decide if to keep new parameters or revert to previous ones
			if (newMSE < currentMSE or math.exp(-params.monteCarloBeta*(newMSE-currentMSE)/currentMSE) > np.random.rand(1)[0]):
				currentMSE, currentRates, currentInitialConditions = newMSE, np.copy(baseRate), np.copy(initialConditions)
			else: # revert to previous parameters
				baseRate, initialConditions = np.copy(currentRates), np.copy(currentInitialConditions)

			# Store the best parameters
			if (newMSE < bestMSE):
				bestMSE, bestRates, bestInitialConditions = newMSE, np.copy(baseRate), np.copy(initialConditions)
				writeBestRateParams(bestMSE, bestRates, bestInitialConditions, outputFileName)
	else: # params.searchMethod == "grid"

		# initializing a vector with nRates + nInitialConditions empty vectors 
		arr = [[] for i in range(params.nRates+params.nInitialConditions)]
		for i in range(0, params.nRates):
			if abs(stepSize[i]) > params.machineEpsilon:
				ratesToTest = np.arange(baseRate[i, 1], baseRate[i, 0], stepSize[i])
			else:
				ratesToTest = np.array(baseRate[i, 1], dtype=float, ndmin=1)
			print(ratesToTest)
			arr[i] = ratesToTest.tolist()
		for i in range(0, params.nInitialConditions):
			if abs(initialConditionStepSize[i]) > params.machineEpsilon:
				initCondsToTest = np.arange(initialConditions[i, 1], initialConditions[i, 0], initialConditionStepSize[i])
			else:
				initCondsToTest = np.array(initialConditions[i, 1], dtype=float, ndmin=1)
			print(initCondsToTest)
			arr[i + params.nRates] = initCondsToTest.tolist()
		
		printCombinationsToFile(arr, "testRatesAndInitialConditions.dat")

		allRatesAndConditions = np.loadtxt("testRatesAndInitialConditions.dat", dtype=float)
		allRates = np.copy(allRatesAndConditions[:,:params.nRates])
		allInitialConditions = np.copy(allRatesAndConditions[:,params.nRates:])
		nIterations = allRatesAndConditions.shape[0]

		# Perform the computationally demanding part

		currentRates = np.zeros((baseRate.shape[0], 1))
		currentInitialConditions = np.zeros((initialConditions.shape[0], 1))
		
		for iIteration in range (0, nIterations):

			currentRates = np.copy(allRates[iIteration])
			currentInitialConditions = np.copy(allInitialConditions[iIteration])

			updateClassicalMasterEquationSystemParFile(currentRates, systemParFileName)

			for j in range (0, params.nInitialConditions):
				
				updateClassicalMasterEquationParamsParFile(currentInitialConditions[j], paramsFileName)
					
				# Copy input files needed for solving the classical master equation
				shutil.copyfile("system.dat", initConditionDirs[j] + "system.par")
				shutil.copyfile("params.dat", initConditionDirs[j] + "params.par")
					
				# Calls c executable kinMCvDt.x 
				allPopulations = solveClassicalMasterEquation(initConditionDirs[j])
					
				# Interpolate data to target data time points and store target population in data
				data[:, j] = np.interp(targetDataTimes, allPopulations[:, 0], allPopulations[:, params.iTargetPopulation])
							
			# Calculate the mean-squared error using scikit-learn library function
			#newMSE = mean_squared_error(normalize(data, axis=0, norm="l2"), normedTargetData)/params.nInitialConditions
			newMSE = mean_squared_error(data, targetData)

			# Write current mean-squared error
			if (iIteration % params.writeIterationStepFrequency == 0):
				writeIteration(iIteration, newMSE, currentRates, currentInitialConditions, iterationFileName)

			# Store the best parameters
			if (newMSE < bestMSE):
				bestMSE, bestRates, bestInitialConditions = newMSE, np.copy(currentRates), np.copy(currentInitialConditions)
				writeBestRateParams(bestMSE, bestRates, bestInitialConditions, outputFileName)

	# Calculate final trajectories using the best fit parameters
	updateClassicalMasterEquationSystemParFile(bestRates, systemParFileName)
	for i in range (0, params.nInitialConditions):
		updateClassicalMasterEquationParamsParFile(bestInitialConditions[i], paramsFileName)
		shutil.copyfile("system.dat", initConditionDirs[i] + "system.par")
		shutil.copyfile("params.dat", initConditionDirs[i] + "params.par")
		allPopulations = solveClassicalMasterEquation(initConditionDirs[i])
		data[:, i] = np.interp(targetDataTimes, allPopulations[:, 0], allPopulations[:, params.iTargetPopulation])

	np.savetxt("bestFitData.dat", np.concatenate((targetDataTimes[:, None], data), axis=1), fmt="%.8f", delimiter=" ")
	shutil.move("bestFitData.dat", outputDir)
	shutil.move("system.dat", outputDir)
	shutil.move("params.dat", outputDir)
	shutil.copy("targetData.par", outputDir + "targetData.dat")
	shutil.move("testRatesAndInitialConditions.dat", outputDir)

###################################################################################

if __name__ == "__main__":
	main()

###################################################################################
