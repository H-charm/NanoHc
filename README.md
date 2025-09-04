# Z peak Validation using CR Z+l

---

## Setup  
Before the first run in each session, load the CMSSW environment and initialize a valid grid proxy:  

```bash
cd CMSSW_13_3_3/src
cmsenv
voms-proxy-init --rfc --voms cms -valid 192:00
```
---

## Running 
To start a run:
```bash
cd PhysicsTools/NanoHc/run  
python3 runHcTrees.py --year <year> --output <output_dir> [--type mc|data] [-n files_per_job]
```
### Arguments
* `--year` : Dataset year (2022, 2022EE, 2023, 2023BPix)
* `--output` : Output directory
* `--type` : "mc" (default) or "data"
* `-n` : Files per job (default = 10, previously 4)
* -`-xsec-file` : Cross-section file (default = samples/xsec.conf)
* `--post` : Merge output files after jobs finish
* `--check-status` : Check submitted jobs progress
* `--check-files` : Check for zombie/incomplete files
* -`-resubmit` : Resubmit failed jobs

---

## Workflow
### 1. Local Testing 
Before submitting to Condor, test changes locally:
```bash
cd jobs_<type>_<year>  
python3 processor.py <job_id>  
```
Note: Check the `metadata.json` file for the job ids 

### 2. Submitting to Condor
Submit jobs with:
```bash
cd jobs_<type>_<year>  
condor_submit submit.sh    
```
For multiple types/years, you can use:
```bash
cd PhysicsTools/NanoHc/run
./submit_all.sh
```
Note: Make sure you have already run ```runHcTrees.py``` for each job before using ```submit_all.sh```.

### 3. Monitoring Jobs

- Check jobs status  
```bash
python3 runHcTrees.py --check-status --year <year> --type <type>
```
- Resubmit failed jobs
```bash
python3 runHcTrees.py --resubmit --year <year> --type <type>
```
- Check for bad output files (zombie or incomplete):
```bash
python3 runHcTrees.py --check-files --year <year> --type <type>
```
Note: You can find log files in each jobs dir   

### 4. Merging Output
After all jobs are finished and validated, merge results:
```bash
python3 runHcTrees.py --post --year <year> --type <type>
```
---

## Modules
You can write new modules in ```python/producers```
Add/remove modules in ```run/static_files/processor.py```
- `jetVetoMapProducer.py`
- `jetJERCProducer.py`
- `electronScaleProducer.py`
- `muonScaleProducer.py`
- `CRProducer.py`
- `puWeightProducer.py`
- `electronSFProducer.py`
- `electronTRGProducer.py`
- `muonSFProducer.py`
- `muonTRGProducer.py`


Add/remove samples in ```run/samples```   


