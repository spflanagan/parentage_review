#!/bin/bash

# Author: Sarah P. Flanagan
# Date: 9 May 2018
# This script runs the parentage simulation program


DATE=`date +%Y%m%d`

###---RUN FROM WHEREVER---###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/parentage_sim"

cd $DIR
cd $PROGDIR

### --- RUN THE PARAMETER COMBINATIONS --- ###
echo "Running the parentage simulation program to generate seven base datasets"
echo "They differ in number of loci, adults, and error rates"
echo "The program will run in the background."
echo "Check the status with htop or by looking at logs/base_parentsim_${DATE}.log"

{
# base 1
./parentsim -L 5000 -S 1 -F 50 -M 50 -o parentsim_F100f4S5000 -d ../../results/

# more adults
./parentsim -L 100 -S 1 -F 250 -M 250 -f 2 -o parentsim_F250f2S100 -d ../../results/
./parentsim -L 100 -S 1 -F 1250 -M 1250 -f 2 -o parentsim_F1250f2S100 -d ../../results/

# allelic dropout
./parentsim -L 100 -S 1 -F 50 -M 50 -f 2 -ad 0.02 -o parentsim_AD2S100 -d ../../results/
./parentsim -L 100 -S 1 -F 50 -M 50 -f 2 -ad 0.05 -o parentsim_AD5S100 -d ../../results/

# scoring errors
./parentsim -L 100 -S 1 -F 50 -M 50 -f 2 -se 0.02 -o parentsim_SE2S100 -d ../../results/
./parentsim -L 100 -S 1 -F 50 -M 50 -f 2 -se 0.05 -o parentsim_SE5S100 -d ../../results/
} >> ../../logs/base_parentsim_${DATE}.log 2>&1 &

