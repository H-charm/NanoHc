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
git clone https://github.com/yiannispar/NanoHc.git PhysicsTools/NanoHc -b dev/nFiles_per_job  
scram b -j8  
```  

Run    
---  
```  
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output dir> -n <files per job> (--post)  
```  

Test locally  
------------  
```  
cd jobs_<year>  
python3 processor.py <job id>  
```  

Submit to condor  
----------------  
```  
cd jobs_<year>  
condor_submit submit.sh  
```  

Merging output files  
--------------------
```  
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output dir> -n <files per job> --post  
```  