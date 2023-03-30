#!/bin/bash

awk 'NR==1{s=$2;next}{printf "%.6f % .8f\n", $1, $2-s; s=$2}' spinVersusBias.dat > dSdVVersusBias.dat
