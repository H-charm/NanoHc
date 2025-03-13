Tree producer for H+c analysis starting from official NANOAOD  

Setup  
-----  
In your **AFS area**:  
```  
cmsrel CMSSW_13_3_3
cd CMSSW_13_3_3/src
cmsenv
export X509_USER_PROXY=/afs/cern.ch/user/${USER:0:1}/$USER/private/x509up_u$(id -u)  
voms-proxy-init --rfc --voms cms -valid 192:00
git clone https://github.com/H-charm/NanoHc.git PhysicsTools/NanoHc
git clone https://github.com/cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
scram b -j8  
```  

Run    
---  
```  
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output dir> (<--type ["mc","data"]>)  
```
Note: It is advised to use the argument ```-n 5``` in order to split the number of files to 5 per job instead of the default 10

Test locally  
------------  
```  
cd jobs_<type>_<year>  
python3 processor.py <job id>  
```
Note: Check the `metadata.json` file for the job ids

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
