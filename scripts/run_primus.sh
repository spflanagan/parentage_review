#!/bin/bash

## Requires already having generated sex and age files in R ##

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../results/"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`

for i in *ped; do
    file=$(echo "$i" | cut -f 1 -d '.')
    ~/Programs/plink --file ${file} --out primus/${file} --genome
    ~/Programs/PRIMUS_v1.9.0/bin/run_PRIMUS.pl --plink_ibd primus/${file}.genome --sex_file primus/${file}.sex --age_file primus/${file}.age
done 2>&1 | tee primus_${DATE}.log

