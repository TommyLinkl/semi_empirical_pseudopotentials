#!/bin/bash

######################################################################################################
# 
# This script creates the input files for a spin-polarized BSE calculation using QE output files
#
# Example use of this script:
# bash ./makeSpinPolarizedParFiles.sh dftOutputFile cubeFileName cubeFilesDirName nHoles nElecs
# bash ./makeSpinPolarizedParFiles.sh temp.out cubeFiles/wfc_K001_B001.cube cubeFiles 12 4
#
# Author: John P. Philbin
# Last modified: November 12th, 2019
#
######################################################################################################


######################################################################################################
# Set useful constants

auToEV=27.2114
evToAU=0.03674930

######################################################################################################
# Grep from $1 (i.e. the quantum espresso output file) the important statistics of the system

nAtoms=( "$(grep 'atoms/cell' $1 | awk '{print $5}')" )
nAtomTypes=( "$(grep 'atomic types' $1 | awk '{print $6}')" )
a=( "$(grep 'lattice parameter (alat)' $1 | awk '{print $5}')" )
cellVolume=( "$(grep 'unit-cell volume' $1 | awk '{print $4}')" )
recipA=( "$(calc 2.0*3.14159265359/$a)" )
nKPoints=( "$(grep 'number of k points=' $1 | awk '{print $5}')" )
nStates=( "$(grep 'number of Kohn-Sham states' $1 | awk '{print $5}')" )
nElectrons=( "$(grep 'number of electrons' $1 | awk '{print $5}')" )
nLines=( "$(awk -v nS=$nStates 'BEGIN {if(nS%8==0){x=1}else{x=2}; print int(nS/8)+x }')" )
homoEnergy=( "$(grep 'highest occupied, lowest unoccupied' $1 | tail -n 1 | awk '{print $7}')" )
lumoEnergy=( "$(grep 'highest occupied, lowest unoccupied' $1 | tail -n 1 | awk '{print $8}')" )
bandGap=( "$(calc $lumoEnergy-$homoEnergy)" )
bandGapAU=( "$(calc $bandGap*$evToAU)" )
fermiEnergy=( "$(calc $homoEnergy+0.5*$bandGap)" )
fermiEnergyAU=( "$(calc $fermiEnergy*$evToAU)" )

######################################################################################################
# Write a configuration file, conf.par

echo $nAtoms > tmpConf1.dat
grep -A"$nAtoms" 'site n.' $1 | tail -n $nAtoms | awk -v l=$a '{printf "%2s % .8f % .8f % .8f\n", $2, $7*l, $8*l, $9*l}' > tmpConf2.dat
cat tmpConf1.dat tmpConf2.dat > conf.par

######################################################################################################
# Determine number of grid points and unit cell used from a .cube ($2) file

nGridPoints=( "$(head -6 $2 | awk '{if(NR>3){print $1}}' | paste -sd" ")" )
nGridX=( "$(echo $nGridPoints | awk '{print $1}')" )
nGridY=( "$(echo $nGridPoints | awk '{print $2}')" )
nGridZ=( "$(echo $nGridPoints | awk '{print $3}')" )
nTotalGridPoints=( "$(head -6 $2 | awk '{if(NR>3){print $1}}' | paste -sd" " | awk '{printf $1*$2*$3}')" )
nPsiLines=( "$(awk -v nTGP=$nTotalGridPoints 'BEGIN {if(nTGP%6==0){x=0}else{x=1}; print int(nTGP/6)+x }')" )
gridPointSpacings=( "$(head -6 $2 | awk '{if(NR>3){print $2" "$3" "$4}}' | paste -sd" ")" )
dX=( "$(calc $a/$nGridX)" )
dY=( "$(calc $a/$nGridY)" )
dZ=( "$(calc $a/$nGridZ)" )
minX=( "$(calc $a*-0.5)" )
minY=( "$(calc $a*-0.5)" )
minZ=( "$(calc $a*-0.5)" )

######################################################################################################
# Write grid.par file 

printf "minPos = %.4f %.4f %.4f\n" "$minX" "$minY" "$minZ" > grid.par

######################################################################################################
# Write a file containing the list of all k-points 

grep ', wk =' $1 | awk -v nKP=$nKPoints '{if(NR<=nKP){printf "% .7f  % .7f  % .7f\n", $5, $6, $7}}' > kpoints_X.par
awk '{printf "%.7f\n", sqrt($1*$1+$2*$2+$3*$3)}' kpoints_X.par > kMagnitudes.dat

######################################################################################################
# Write a file containing the band energies for all the k-points 

grep -A"$nLines" 'bands (ev):' $1 > tmp1.dat
awk -v nL=$nLines '{if($1=="k"){x=0;y=NR;}else{x++;y++;};if((x>1)&&(x<=nL)){print $0} }' tmp1.dat > tmp2.dat
nLines=( "$(awk -v nS=$nStates 'BEGIN {if(nS%8==0){x=0}else{x=1}; print int(nS/8)+x }')" )
awk -v nL=$nLines '{line=line " " $0} (NR%nL)==0{print substr(line,2); line=""}' tmp2.dat > tmp3.dat 
awk -v nS=$nStates '{for(i=1;i<=nS;i++){printf "% .12f ", $i; if(i==nS){printf "\n"}} }' tmp3.dat > bandEnergies.dat
paste -d" " kMagnitudes.dat bandEnergies.dat > expBandStruct_X.par

######################################################################################################
# Write a file that contains an ordered list of all the energies, irrespective of the k-points

awk -v nS=$nStates '{for(i=1;i<=nS;i++){printf "% .12f\n", $i} }' bandEnergies.dat > tmp4.dat 
sort -n tmp4.dat > tmp5.dat
awk '{printf "%5d % .12f\n", x, $1; x++ }' tmp5.dat > allEnergies.par
awk -v nB=$nStates '{if(NR<=nB){print $0}}' tmp4.dat > tmp6.dat 
awk '{printf "%5d   0.5  % .12f\n", x, $1; x++ }' tmp6.dat > spinUpAllEnergies.par
awk -v nB=$nStates '{if(NR>nB){print $0}}' tmp4.dat > tmp6.dat
awk '{printf "%5d  -0.5  % .12f\n", x, $1; x++ }' tmp6.dat > spinDownAllEnergies.par
cat spinUpAllEnergies.par spinDownAllEnergies.par | sort -nk 3 > tmp7.dat
awk '{printf "%5d % .1f % .12f\n", x, $2, $3; x++ }' tmp7.dat > allEnergiesWithSpin.par
awk '{printf "%5d % .1f % .12f\n", $1+1, $2, $3;}' tmp7.dat > allEnergiesWithSpinPsiIndices.par
magnetization=( "$(awk -v eF=$fermiEnergy '{if($3<eF){x+=$2}}END{print x}' allEnergiesWithSpin.par)" )

######################################################################################################
# Write allEval.par file -> add very low sigma values

awk -v c=$evToAU -v s=0.00000001 '{printf "%5d % .12f  %.12f\n", $1, $2*c, s}' allEnergies.par > allEval.par
awk -v c=$evToAU -v s=0.00000001 '{printf "%5d % .1f % .12f  %.12f\n", $1, $2, $3*c, s}' allEnergiesWithSpin.par > allSpinEval.par
awk -v c=$evToAU -v s=0.00000001 '{printf "%5d % .1f % .12f  %.12f\n", $1, $2, $3*c, s}' allEnergiesWithSpinPsiIndices.par > allSpinEvalPsiIndices.par

######################################################################################################
# Determine indices of HOMO, LUMO, highest energy hole and highest energy electron 

iHOMO=( "$(awk -v nE=$nElectrons 'BEGIN {print int(nE/2 + nE%2)}')" )
iMinHole=( "$(awk -v nE=$nElectrons -v nPsiH=$4 'BEGIN {print int(nE/2 + nE%2 + 1.0 - nPsiH)}')" )
iLUMO=( "$(awk -v nE=$nElectrons 'BEGIN {print int(nE/2 + nE%2 + 1.0)}')" )
iMaxElec=( "$(awk -v nE=$nElectrons -v nPsiE=$5 'BEGIN {print int(nE/2 + nE%2 + nPsiE)}')" )
nFinalStates=( "$(calc $iMaxElec-$iMinHole+1)")
iMinHoleEval=( "$(calc $nElectrons-$4)")
iMaxElecEval=( "$(calc $nElectrons+$5-1)")

######################################################################################################
# Make final eval.par and spinEval.par files

awk -v h=$iMinHoleEval -v e=$iMaxElecEval '{if($1 >= h && $1 <= e){print $0}}' allEval.par > tmpEval.dat
awk '{$1=(NR-1); printf "%4d % .12f  %.12f\n", $1, $2, $3;}' tmpEval.dat > 'eval'.par
awk -v h=$iMinHoleEval -v e=$iMaxElecEval '{if($1 >= h && $1 <= e){print $0}}' allSpinEval.par > tmpSpinEval.dat
awk '{$1=(NR-1); printf "%4d % .1f % .12f  %.12f\n", $1, $2, $3, $4;}' tmpSpinEval.dat > spinEval.par

######################################################################################################
# Order the names of the wavefunction cube files based on order and indices 

awk -v h=$iMinHoleEval -v e=$iMaxElecEval -v d=$3 '{if(NR-1 >= h && NR-1 <= e){ if($2>0.01){x="1"}else{x="2"}; printf "%s/wfc_K00%1d_B%04d.cube\n", d, x, $1; } }' allSpinEvalPsiIndices.par > psiFilenames.par

######################################################################################################
# Delete files that will be concatenated

if [ -e "psi.par" ]; then
	rm psi.par
fi
if [ -e "infoPsi.par" ]; then
	rm infoPsi.par
fi

######################################################################################################
# Write the final psi.par file

for i in $(seq 1 $nFinalStates); do
	cubeFileName=( "$(awk -v n=$i '{if(NR==n){print $0}}' psiFilenames.par)" )
	echo $cubeFileName >> infoPsi.par
	tail -n $nPsiLines $cubeFileName | fmt -1 > tmpPsi.dat
	awk '{printf "%.16f\n", $1}' tmpPsi.dat > tmpPsi2.dat
	cat tmpPsi2.dat >> psi.par
done

######################################################################################################
# Delete temporary files

rm tmp*.dat kMagnitudes.dat bandEnergies.dat

######################################################################################################
# Write important statistics from the calculation

printf "Number of atoms       = %d\n"    "$nAtoms"
printf "Number of atom types  = %d\n"    "$nAtomTypes"
printf "Number of electrons   = %.1f\n"  "$nElectrons"
printf "Number of grid points = %d\n"    "$nTotalGridPoints"
printf "Grid points           = %s\n"    "$nGridPoints"
printf "gridPointSpacings     = %.6f %.6f %.6f\n"  "$dX" "$dY" "$dZ"
printf "Lattice parameter     = %.6f\n"  "$a"
printf "Box min grid point    = %.6f %.6f %.6f\n"  "$minX" "$minY" "$minZ"
printf "Unit-cell volume      = %.2f\n"  "$cellVolume"
printf "k-space scale         = %.6f\n"  "$recipA"
printf "Number of k-points    = %d\n"    "$nKPoints"
printf "Number of bands       = %d\n"    "$nStates"
printf "HOMO index            = %d\n"    "$iHOMO"
printf "LUMO index            = %d\n"    "$iLUMO"
printf "HOMO energy           = %.6f\n"  "$homoEnergy"
printf "LUMO energy           = %.6f\n"  "$lumoEnergy"
printf "Band gap              = %.6f\n"  "$bandGap"
printf "Band gap              = %.6f\n"  "$bandGapAU"
printf "Fermi energy          = %.6f\n"  "$fermiEnergy"
printf "Fermi energy          = %.8f\n"  "$fermiEnergyAU" 
printf "Minumum hole index    = %d\n"    "$iMinHole"
printf "Maximum elec index    = %d\n"    "$iMaxElec"
printf "Magnetization         = %.3f\n"  "$magnetization"
printf "Number of quasiholes  = %d\n"    "$4"
printf "Number of quasielecs  = %d\n"    "$5"

######################################################################################################