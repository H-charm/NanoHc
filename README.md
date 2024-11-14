Tree producer for H+c analysis starting from official NANOAOD  

Setup  
-----  
In your **AFS area**:  
```  
cmsrel CMSSW_13_3_0  
cd CMSSW_13_3_0/src  
cmsenv  
export X509_USER_PROXY=/afs/cern.ch/user/${USER:0:1}/$USER/private/x509up_u$(id -u)  
voms-proxy-init --rfc --voms cms -valid 192:00  
git clone https://github.com/H-charm/NanoHc.git PhysicsTools/NanoHc   
scram b -j8  
```  

Run    
---  
```  
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output dir> (<--type ["mc","data"]>)  
```  

Test locally  
------------  
```  
cd jobs_<type>_<year>  
python3 processor.py <job id>  
```  

Submit to condor  
----------------  
```  
cd jobs_<type>_<year>  
condor_submit submit.sh    
```  

Check jobs status  
----------------  
Run again ```runHcTrees.py``` with ```--check-status```  

Merging output files  
--------------------
Run again ```runHcTrees.py``` with ```--post```  

Resubmit failed jobs    
--------------------  
Run again ```runHcTrees.py``` with ```--resubmit```  
Go to jobs dir and submit again to condor  

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
- ```--resubmit```: Resubmit failed jobs  