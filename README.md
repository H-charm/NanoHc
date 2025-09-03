# Z peak Validation

## Run

Before running for the first time (for each session) you need to load the CMSSW for the required packages and to have a valid grid proxy for the NanoAOD samples.

```bash
cd CMSSW_13_3_3/src
cmsenv
voms-proxy-init --rfc --voms cms -valid 192:00
```
In order to run this repository you need to do:
```
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output dir> (<--type ["mc","data"]>)  
```
In addition you can you use the argument ```-n X``` in order to split the number of files to X per job instead of the default 4

### Test locally  
If any changes are being made, it is advised to test them locally before submitting to condor. This can be done by
```  
cd jobs_<type>_<year>  
python3 processor.py <job id>  
```
Note: Check the `metadata.json` file for the job ids 

### Submit to condor  
In order to submit the jobs:
```  
cd jobs_<type>_<year>  
condor_submit submit.sh    
```  
Note: You can use the `run/submit_all.sh` script in order to submit multiple types,years requiring that you have run the `runHcTrees.py` script before for every job.


Check jobs status  
----------------  
Run again ```runHcTrees.py``` with ```--check-status``` in order to check the progress of the submitted jobs.

If some jobs have failed for some reason, you can resubmit them by running again ```runHcTrees.py``` with ```--resubmit```.

Check file status
----------------  
Run again ```runHcTrees.py``` with ```--check-files``` in order to check for Zombie or incomplete files.

Merging output files  
--------------------
Run again ```runHcTrees.py``` with ```--post``` in order to merge.

## Producers

Important notes 
--------------  
- Add/remove modules in ```run/static_files/processor.py```  
- Add/remove samples in ```run/samples```   
- You can find log files in each jobs dir   
- You can write new modules in ```python/producers```  

Argument full list  
------------------  
- ```--year```    
- ```--output```    
- ```--type```: "mc" or "data", default = "mc"  
- ```--post```: Merge output files  
- ```-n```: Number of files per job, default = 10 
- ```--xsec-file```: xsec file, default = "samples/xsec.conf"  
- ```--check-status```: Check jobs status  

