#!/bin/bash

######################################################################################################
# 
# This script creates the input files for a spin-polarized BSE calculation using QE output cube files
# and an energy file containing the eigenenergies and occupation numbers
#
# Example use of this script:
# bash ./makeSpinPolarizedParFilesGWEnergies.sh energiesFile cubeFilesDirName nHoles nElecs qeDFTOutputFileName
# bash ./makeSpinPolarizedParFilesGWEnergies.sh GW_energies.dat cubeFiles 12 4 temp.out
#
# Author: John P. Philbin
# Last modified: November 14th, 2019
#
######################################################################################################


######################################################################################################
# Set useful constants

auToEV=27.2114
evToAU=0.03674930

######################################################################################################
# Make 

#sort -nk 3 $1 | awk '{printf "%d % .3f % .1f % .8f\n", NR, $1,  0.5, $3}' > dftSpinUpEnergies.dat
#sort -nk 4 $1 | awk '{printf "%d % .3f % .1f % .8f\n", NR, $2, -0.5, $4}' > dftSpinDownEnergies.dat
#cat dftSpinDownEnergies.dat dftSpinUpEnergies.dat | sort -nk 4 | awk '{printf "%ld %s\n", NR, $0}'> allDFTEnergies.dat

sort -nk 5 $1 | awk '{printf "%d % .3f % .1f % .8f\n", NR, $1,  0.5, $5}' > gwSpinUpEnergies.dat
sort -nk 6 $1 | awk '{printf "%d % .3f % .1f % .8f\n", NR, $2, -0.5, $6}' > gwSpinDownEnergies.dat
cat gwSpinDownEnergies.dat gwSpinUpEnergies.dat | sort -nk 4 | awk '{printf "%ld %s\n", NR, $0}'> allGWEnergies.dat

# Identify HOMO and LUMO
iHOMO=( "$(awk '{if($3>0.5){x=$0}}END{print x}' allGWEnergies.dat | awk '{print $1}' )")
homoEnergy=( "$(awk '{if($3>0.5){x=$0}}END{print x}' allGWEnergies.dat | awk '{print $5}' )")
iLUMO=( "$(awk '{if($3<0.5){print $0; exit }}' allGWEnergies.dat | awk '{print $1}' )")
lumoEnergy=( "$(awk '{if($3<0.5){print $0; exit }}' allGWEnergies.dat | awk '{print $5}' )")
bandGap=( "$(calc $lumoEnergy-$homoEnergy)" )
bandGapAU=( "$(calc $bandGap*$evToAU)" )
fermiEnergy=( "$(calc $homoEnergy+0.5*$bandGap)" )
fermiEnergyAU=( "$(calc $fermiEnergy*$evToAU)" )

# Determine the number of electrons
nElectrons=( "$(awk '{if($3>0.5){x+=1}}END{print x}' allGWEnergies.dat )")

# Determine the total magnetization
totalMagnetization=( "$(awk '{if($3>0.5){x+=$4}}END{print x}' allGWEnergies.dat )")

######################################################################################################
# Make final eval.par and spinEval.par files

awk -v nH=$3 -v nE=$4 -v iH=$iHOMO -v iL=$iLUMO '{if( (NR>(iH-nH)) && (NR<(iL+nE)) ){print $0}}' allGWEnergies.dat > selectedGWEnergies.dat
awk '{printf "%5d % .1f % .12f % .8f\n", NR-1, $4, $5, 0.00000010}' selectedGWEnergies.dat > tmp2.dat
awk -v s=$evToAU '{printf "%5d % .1f % .12f % .8f\n", NR-1, $4, $5*s, 0.00000010}' selectedGWEnergies.dat > spinEval.par
awk '{printf "%5d % .12f % .8f\n", $1, $3, $4}' spinEval.par > eval.par

# Determine the magnetization
magnetization=( "$(awk -v eF=$fermiEnergyAU '{if($3<eF){x+=$2}}END{print x}' spinEval.par )")

######################################################################################################
# Determine number of grid points and unit cell used from a .cube file

# Print a file containing the list of cube files from which the psi.par will be built
awk -v dir=$2 '{if($4>0.0){x="1"}else{x="2"}; printf "%s/wfc_K00%1d_B%04d.cube\n", dir, x, $2;}' selectedGWEnergies.dat > psiFilenames.dat

cubeFileName=( "$(awk '{if(NR==1){print $0}}' psiFilenames.dat)" )
nTotalGridPoints=( "$(head -6 $cubeFileName | awk '{if(NR>3){print $1}}' | paste -sd" " | awk '{printf $1*$2*$3}')" )
nAtoms=( "$(awk '{if(NR==3){print $1; exit}}' $cubeFileName)" )
nTotalLines=( "$(wc -l $cubeFileName | awk '{print $1}')" )
nPsiLines=( "$(calc $nTotalLines-$nAtoms-6)" )

######################################################################################################
# Write a configuration file, conf.par

echo $nAtoms > tmpConf1.dat
a=( "$(grep 'lattice parameter (alat)' $5 | awk '{print $5}')" )
grep -A"$nAtoms" 'site n.' $5 | tail -n $nAtoms | awk -v l=$a '{printf "%2s % .8f % .8f % .8f\n", $2, $7*l, $8*l, $9*l}' > tmpConf2.dat
cat tmpConf1.dat tmpConf2.dat > conf.par

######################################################################################################
# Make grid.par file 

nGridPoints=( "$(head -6 $cubeFileName | awk '{if(NR>3){print $1}}' | paste -sd" ")" )
nGridX=( "$(echo $nGridPoints | awk '{print $1}')" )
nGridY=( "$(echo $nGridPoints | awk '{print $2}')" )
nGridZ=( "$(echo $nGridPoints | awk '{print $3}')" )
gridPointSpacings=( "$(head -6 $cubeFileName | awk '{if(NR>3){print $2" "$3" "$4}}' | paste -sd" ")" )
dX=( "$(awk '{if(NR==4){print $2; exit }}' $cubeFileName)" )
dY=( "$(awk '{if(NR==5){print $3; exit }}' $cubeFileName)" )
dZ=( "$(awk '{if(NR==6){print $4; exit }}' $cubeFileName)" )
minX=( "$(calc $dX*$nGridX*-0.5)" )
minY=( "$(calc $dY*$nGridY*-0.5)" )
minZ=( "$(calc $dZ*$nGridZ*-0.5)" )

# Write grid.par file 
printf "minPos = %.4f %.4f %.4f\n" "$minX" "$minY" "$minZ" > grid.par

######################################################################################################
# Write the final psi.par file

# Delete file that will be concatenated
if [ -e "psi.par" ]; then
	rm psi.par
fi

nFinalStates=( "$(calc $3+$4)")
for i in $(seq 1 $nFinalStates); do
	cubeFileName=( "$(awk -v n=$i '{if(NR==n){print $0}}' psiFilenames.dat)" )
	tail -n $nPsiLines $cubeFileName | fmt -1 > tmpPsi.dat
	awk '{printf "%.16f\n", $1}' tmpPsi.dat > tmpPsi2.dat
	cat tmpPsi2.dat >> psi.par
done

######################################################################################################
# Write important statistics from the calculation

printf "Number of atoms          = %d\n"    "$nAtoms"
printf "Number of electrons      = %d\n"    "$nElectrons"
printf "Number of grid points    = %d\n"    "$nTotalGridPoints"
printf "Grid points              = %s\n"    "$nGridPoints"
printf "gridPointSpacings        = %.6f %.6f %.6f\n"  "$dX" "$dY" "$dZ"
printf "Box min grid point       = %.6f %.6f %.6f\n"  "$minX" "$minY" "$minZ"
printf "HOMO index               = %d\n"    "$iHOMO"
printf "LUMO index               = %d\n"    "$iLUMO"
printf "HOMO energy              = %.8f\n"  "$homoEnergy"
printf "LUMO energy              = %.8f\n"  "$lumoEnergy"
printf "Band gap                 = %.8f\n"  "$bandGap"
printf "Band gap                 = %.10f\n" "$bandGapAU"
printf "Fermi energy             = %.8f\n"  "$fermiEnergy"
printf "Fermi energy             = %.10f\n" "$fermiEnergyAU" 
printf "System magnetization     = %.3f\n"  "$totalMagnetization"
printf "Subsystem magnetization  = %.3f\n"  "$magnetization"
printf "Number of quasiholes     = %d\n"    "$3"
printf "Number of quasielecs     = %d\n"    "$4"
printf "Number of quasiparticles = %d\n"    "$nFinalStates"

######################################################################################################
# Delete temporary files

rm tmp*.dat 

######################################################################################################