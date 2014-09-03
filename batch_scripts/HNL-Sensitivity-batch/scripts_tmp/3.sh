#!/bin/bash

export RFIO_USE_CASTOR_V2=YES
export STAGE_HOST=castorpublic
export PYTHIA8=/afs/cern.ch/user/e/egraveri/pythia8186
export PYTHIA8DATA=$PYTHIA8/xmldoc
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
source /afs/cern.ch/user/e/egraveri/my-root/root_5.34.20/bin/thisroot.sh
export PYTHONDIR=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc48-opt
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH

cd /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/

python /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/sandbox/scan.py  3 /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/sandbox  > /afs/cern.ch/work/e/egraveri/public/SHiP/HNL-Sensitivity-batch/logs/3-6.log 
