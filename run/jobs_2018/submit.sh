
executable = condor_exec.sh

arguments = $(jobid) 

request_memory  = 2000
request_disk    = 10000000

error   = log/err.$(Process)
output  = log/out.$(Process)
log     = log/logFile.log

+JobFlavour = "longlunch"

queue jobid from job_ids.txt
    