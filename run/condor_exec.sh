#! /usr/bin/env sh

## arguments
jobid=$1
cmssw_base=$2
jobs_dir=$3

## needed to load CMSSW libraries/packages
cd $cmssw_base/src
eval `scramv1 runtime -sh`

export X509_USER_PROXY=/afs/cern.ch/user/i/iparaske/private/x509up_u147242

## run executable
cd $jobs_dir
python3 processor.py $jobid