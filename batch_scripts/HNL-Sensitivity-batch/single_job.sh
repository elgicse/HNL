#!/bin/bash

CURRENTMODEL=3

cd /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch
cp base.sh scripts_tmp/${CURRENTMODEL}.sh
chmod +x scripts_tmp/${CURRENTMODEL}.sh
echo "python /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/sandbox/scan.py "  $CURRENTMODEL $PWD/sandbox  " > /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/logs/${CURRENTMODEL}.log " >> scripts_tmp/${CURRENTMODEL}.sh
bsub -q test -W 0:20 -J eg-SHIP-${CURRENTMODEL} /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/scripts_tmp/${CURRENTMODEL}.sh