#!/bin/bash

MODELS=3
CURRENTMODEL=1

for ((i=0 ; i<$MODELS; ++i ))
do

    cd /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch
    cp base.sh scripts_tmp/${CURRENTMODEL}.sh
    chmod +x scripts_tmp/${CURRENTMODEL}.sh
    echo "python /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/sandbox/scan.py " $CURRENTMODEL $PWD/sandbox  " > /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/logs/${CURRENTMODEL}-6.log " >> scripts_tmp/${CURRENTMODEL}.sh
    # ./scripts_tmp/${i}.sh > logs/{$i}_{$SEED}.log
#    bsub -q 8nm -o /dev/null -e /dev/null -J eg-SHIP-${CURRENTMODEL} /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/scripts_tmp/${CURRENTMODEL}.sh
    bsub -q 2nw -J eg-SHIP-${CURRENTMODEL} /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/scripts_tmp/${CURRENTMODEL}.sh
    CURRENTMODEL=$(( CURRENTMODEL + 1 ))
      
done