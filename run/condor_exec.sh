#! /usr/bin/env sh

## arguments
jobid=$1

## needed to load CMSSW libraries/packages
cd /afs/cern.ch/user/i/iparaske/test/CMSSW_13_3_0/src/PhysicsTools/NanoHc/run/
eval `scramv1 runtime -sh`

export X509_USER_PROXY=/afs/cern.ch/user/i/iparaske/private/x509up_u147242

## run executable
cd /afs/cern.ch/user/i/iparaske/test/CMSSW_13_3_0/src/PhysicsTools/NanoHc/run/jobs_2018/
python3 processor.py $jobid