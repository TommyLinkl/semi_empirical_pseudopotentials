#!/bin/bash

head -2 projRhoHole-ni-0.dat | tail -n 1 | awk '{printf "#BoxSize % .3f % .3f % .3f\n#HalfBox % .3f % .3f % .3f\n#i  xAverage yAverage zAverage\n", -$1, -$3, -$5, -$1/2.0, -$3/2.0, -$5/2.0}'

for i in {0..39}
do
	awk -v c=$i '{if(NR>1){x+=($1*$2);y+=($3*$4);z+=($5*$6)}}END{printf "%3d % 8.3f % 8.3f % 8.3f\n", c, x, y, z}' projRhoHole-ni-$i.dat
done
