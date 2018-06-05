#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../results/SNPPIT/"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`

for i in *snppit.txt; do
    echo ${i}
    ./snppit -f ${i}
    for f in snppit_output*; do mv "$f" "${i%snppit.txt}${f}"; done
done 2>&1 | tee snppit_${DATE}.log

