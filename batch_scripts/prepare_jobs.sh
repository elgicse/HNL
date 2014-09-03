#!/bin/bash

NEVENTS=30000000
SEED=5000
NJOBS=500

for ((i=0 ; i<$NJOBS; ++i ))
do

    cd /afs/cern.ch/work/e/egraveri/public/SHiP/batch_scripts    
    cp base.sh scripts_tmp/${SEED}.sh
    chmod +x scripts_tmp/${SEED}.sh
    echo "python2.6 /afs/cern.ch/work/e/egraveri/public/SHiP/batch_scripts/pythia8_Beauty.py"  $NEVENTS $SEED  " > /afs/cern.ch/work/e/egraveri/public/SHiP/batch_scripts/logs/${SEED}.log " >> scripts_tmp/${SEED}.sh
    # ./scripts_tmp/${i}.sh > logs/{$i}_{$SEED}.log
    bsub -q 2nd -o /dev/null -e /dev/null -J eg-b-${SEED} /afs/cern.ch/work/e/egraveri/public/SHiP/batch_scripts/scripts_tmp/${SEED}.sh
    SEED=$(( SEED + 1 ))
      
done
