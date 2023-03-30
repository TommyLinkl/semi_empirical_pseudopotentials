#!/bin/bash
#
# $1 = csv file name
# $2 = optical gap (final electronic energy for all systems (E(t>1ns) = Eopt))
# $3 = pump energy (energy of one initial exciton created)
# $4, $5, ... = number of initial excitons created by the pump
#

let "nPumpIntensites = $# - 3"

# Print csv with space delimeters and no header line
awk -F, '{if(NR>1){printf "% .1f ", $1; for (i = 2; i <=NF; i++) printf  "% .5f ", $i; printf "\n"}}' $1 > tmpCSVToTargetData.par

# Make sure the lowest times are on top
sort -nk 1 tmpCSVToTargetData.par > tmpCSVToTargetData2.par

# Remove times before time equals 0
awk '{if($1>-0.001){print $0}}' tmpCSVToTargetData2.par > tmpCSVToTargetData3.par

# Average the last few times in order to determine the level at which each power intensity plateaus.
awk '{if($1>949.){n+=1.0; for (i = 2; i <= NF; i++) x[i-2]+=$i}}END{for (i = 0; i <= NF-2; i++){ print x[i]/n" "1.0-x[i]/n; a+=x[i]}; print n" "NF-1" "a/n/(NF-1.)}' tmpCSVToTargetData3.par > tmpCSVToTargetData4.par 

c=0.0
for (( i=4; i<=$#; i++ )); do
	t=( "$(awk -v l=$i -v n=${!i} -v pE=$3 -v eF=$2 '{if(NR==(l-3)){x=$2/(pE*n-eF)}}END{print x}' tmpCSVToTargetData4.par)" )
	c=$( python -c "print $c + $t" )
	let "index = $i - 4"
	tmp=$( python -c "print 1.0*${!i}" )
	array[$index]=$tmp
	cArray[$index]=$t
done
c=$( python -c "print $c/$nPumpIntensites" )

#c=0.00425780666667

awk -v s=$c -v pE=$3 -v n0="${array[*]}" 'BEGIN{split(n0, a, " ")}{if($1>-0.001){printf "% .1f ", $1; for (i = 2; i <= NF; i++) printf "% .8f ", pE*a[i-1]-((1.0-$i)/s); printf "\n"}}' tmpCSVToTargetData3.par > targetData.par

awk -v s="${cArray[*]}" -v pE=$3 -v n0="${array[*]}" 'BEGIN{split(n0, a, " "); split(s, b, " ")}{if($1>-0.001){printf "% .1f ", $1; for (i = 2; i <= NF; i++) printf "% .8f ", pE*a[i-1]-((1.0-$i)/b[i-1]); printf "\n"}}' tmpCSVToTargetData3.par >  targetData_manyCs.par

echo csv file = $1
echo number of pump intensities = $nPumpIntensites
echo optical gap = $2
echo pump energy = $3
for (( i=4; i<=$#; i++ )); do
	echo "number of initial excitons = ${!i}"
done
echo ${cArray[*]}
echo cAverage = $c

#rm tmpCSVToTargetData*.par
