#!/bin/bash

COUNTER=60                       #< count the current progression
TOTAL=100
#TOTAL=` ${indir} -name "*.zip" | wc -l`  #< total number of files
while [ ${COUNTER} -lt ${TOTAL} ]
do
    echo "Processing file ${COUNTER} / ${TOTAL}"
   
    _JOB=`python [3]healpix_count.py ${COUNTER}`

    echo ${_JOB}
    #echo ${COUNTER}
    COUNTER=$(( COUNTER + 1 ));   
done
