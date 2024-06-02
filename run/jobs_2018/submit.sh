
executable = condor_exec.sh

arguments = $(jobid) 

error   = log/err.$(Process)
output  = log/out.$(Process)
log     = log/logFile.log

+JobFlavour = "longlunch"

queue jobid from job_ids.txt
    