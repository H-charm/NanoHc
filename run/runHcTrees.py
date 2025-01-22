#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import yaml
import json
import os
import sys
from pathlib import Path
import helpers

## parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--year', type=str, help='Year to run', default="2022")
parser.add_argument('--output', type=str, help='Output dir', default = "/eos/user/i/iparaske/HcTrees/")
parser.add_argument('--type', type=str, help='mc or data', default = "mc", choices=['mc', 'data'])
parser.add_argument('--post',help='Merge output files',action='store_true')
parser.add_argument('-n',type=int, help='Number of files per job', default=10)
parser.add_argument('--xsec-file', type=str, help='xsec file', default = "samples/xsec.conf")
parser.add_argument('--check-status', help='Checks jobs status', action='store_true')
parser.add_argument('--resubmit', help='Resubmit failed jobs', action='store_true')
args = parser.parse_args()

golden_json = {
    '2022': 'Cert_Collisions2022_355100_362760_Golden.json',
}

## read samples yaml file and produce json file to be used by condor
def create_metadata_json():
    
    dataset_type = args.type
    jobs_dir = "jobs_" + dataset_type + "_" + args.year
    year = args.year
    
    ## read samples yaml file
    samples_yaml_file = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoHc/run/samples/" + dataset_type + "_" + year + ".yaml"
    with open(samples_yaml_file, 'r') as file:
        samples = yaml.safe_load(file)

    physics_processes = []
    das_dict = {}
    for sample in samples:        
        das_dict[sample] = {}
        for dataset in samples[sample]:
            ## find files using DAS
            print("DAS query for dataset " + dataset)
            das_query = 'dasgoclient --query="file dataset=' + dataset + '"'
            query_out = os.popen(das_query)
            files_found = ['root://xrootd-cms.infn.it/'+_file.strip() for _file in query_out]
            physics_process = dataset.split("/")[1]
            physics_processes.append(physics_process)
            das_dict[sample][physics_process] = files_found
            print(f"{len(files_found)} files found")
                                    
    ## write json file
    json_file = jobs_dir + '/metadata.json'
    json_content = {}
    
    json_content["output_dir"] = args.output
    json_content["jobs_dir"] = args.jobs_dir
    json_content["year"] = args.year    
    json_content["type"] = dataset_type
    if dataset_type == "data":  
        json_content["golden_json"] = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoHc/data/JSON/" + golden_json[args.year]
    else:
        json_content["golden_json"] = None
    json_content["sample_names"] = []
    json_content["physics_processes"] = []
    json_content["jobs"] = []
        
    for sample_name in samples: json_content["sample_names"].append(sample_name)
    for physics_process in physics_processes: json_content["physics_processes"].append(physics_process)
    
    job_id = 0
    for sample in samples:
        for physics_process in das_dict[sample]:
            for chunk in enumerate(helpers.get_chunks(das_dict[sample][physics_process],args.n)):
                json_content["jobs"].append({"job_id": job_id, "input_files": chunk[1], "sample_name": sample ,"physics_process": physics_process})
                job_id += 1

    with open(json_file, 'w') as file:
        json.dump(json_content, file, indent=4)
        
def write_condor_submit(jobids_file):
    
    cmssw_base = os.environ['CMSSW_BASE']
    jobs_dir_path =os.getcwd() + "/" + args.jobs_dir 
    
    condor_submit_file = open(args.jobs_dir + "/submit.sh","w")
    condor_submit_file.write('''
executable = condor_exec.sh

arguments = $(jobid) ''' + cmssw_base + ''' ''' + jobs_dir_path  + ''' 

request_memory  = 2000
request_disk    = 10000000

output = log/$(jobid).out
error = log/$(jobid).err
log = log/$(jobid).log

JobBatchName = HcTrees_''' + args.type + '''_''' + args.year + '''
+JobFlavour = "tomorrow"

queue jobid from ''' + jobids_file)

    condor_submit_file.close()

def create_condor_submit():
    
    ## get job ids from json
    with open(args.jobs_dir + "/metadata.json", 'r') as file:
        data = json.load(file)  
    njobs = len(data["jobs"])
    jobs_list = [i for i in range(njobs)]
    
    ## write job ids to txt file
    with open(args.jobs_dir + "/job_ids.txt", 'w') as file:
        for index, item in enumerate(jobs_list):
            if index < len(jobs_list) - 1:
                file.write(str(item) + '\n')
            else:
                file.write(str(item))
    
    write_condor_submit(jobids_file="job_ids.txt")

def merge_output_files():

    file_path = args.jobs_dir + '/metadata.json'
    with open(file_path, 'r') as file:
        data = json.load(file)
        
    base_output_dir = data["output_dir"]
    dataset_type = data["type"]
    year = data["year"]
    merged_dir = os.path.join(base_output_dir, dataset_type, year, "merged")
    os.system("mkdir -p " + merged_dir)
    
    sample_dirs = [d.name for d in Path(os.path.join(base_output_dir, dataset_type, year)).iterdir() if d.is_dir()]
    
    for sample_dir in sample_dirs:
        if sample_dir == "merged": continue
        input_files = []
        output_file = os.path.join(merged_dir, sample_dir + "_tree.root")
        
        physics_process_dirs = [d.name for d in Path(os.path.join(base_output_dir, dataset_type, year, sample_dir)).iterdir() if d.is_dir()]

        for physics_process_dir in physics_process_dirs:
            infiles_dir = os.path.join(base_output_dir, dataset_type, year, sample_dir, physics_process_dir)
            
            for root, _, files in os.walk(infiles_dir):
                for file in files:
                    input_files.append(os.path.join(root, file))
                          
        cmd = 'haddnano.py {outfile} {infiles}'.format(outfile=output_file, infiles=' '.join(input_files))
        os.system(cmd)

def parse_sample_xsec(cfgfile):
    xsec_dict = {}
    with open(cfgfile) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            pieces = l.split()
            samp = None
            xsec = None
            isData = False
            for s in pieces:
                if '/MINIAOD' in s or '/NANOAOD' in s:
                    samp = s.split('/')[1]
                    if '/MINIAODSIM' not in s and '/NANOAODSIM' not in s:
                        isData = True
                        break
                else:
                    try:
                        xsec = float(s)
                    except ValueError:
                        try:
                            import numexpr
                            xsec = numexpr.evaluate(s).item()
                        except:
                            pass
            if samp is None:
                print('Ignore line:\n%s' % l)
            elif not isData and xsec is None:
                print('Cannot find cross section:\n%s' % l)
            else:
                if samp in xsec_dict and xsec_dict[samp] != xsec:
                    raise RuntimeError('Inconsistent entries for sample %s' % samp)
                xsec_dict[samp] = xsec
                if 'PSweights_' in samp:
                    xsec_dict[samp.replace('PSweights_', '')] = xsec
                    
    return xsec_dict

def add_weights(file, xsec, lumi=1000., treename='Events'):
    from array import array
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    def _get_sum(tree, wgtvar):
        htmp = ROOT.TH1D('htmp', 'htmp', 1, 0, 10)
        tree.Project('htmp', '1.0', wgtvar)
        return float(htmp.Integral())
    
    def _fill_const_branch(tree, branch_name, buff, lenVar=0):
        if lenVar > 0:
            b = tree.Branch(branch_name, buff, branch_name + "[" + str(nScaleWeights) + "]" + '/F')
            b.SetBasketSize(tree.GetEntries() * 2)  # be sure we do not trigger flushing
            for _ in range(tree.GetEntries()):
                b.Fill()
        else:
            b = tree.Branch(branch_name, buff, branch_name + '/F')
            b.SetBasketSize(tree.GetEntries() * 2)  # be sure we do not trigger flushing
            for _ in range(tree.GetEntries()):
                b.Fill()
                
    f = ROOT.TFile(str(file), 'UPDATE')
    run_tree = f.Get('Runs')
    tree = f.Get(treename)

    # fill cross section weights to the 'Events' tree
    sumwgts = _get_sum(run_tree, 'genEventSumw')
    xsecwgt = xsec * lumi / sumwgts
    xsec_buff = array('f', [xsecwgt])
    _fill_const_branch(tree, "xsecWeight", xsec_buff)

    # fill LHE weight re-normalization factors
    if tree.GetBranch('LHEScaleWeight'):
        run_tree.GetEntry(0)
        nScaleWeights = run_tree.nLHEScaleSumw
        scale_weight_norm_buff = array('f',
                                       [sumwgts / _get_sum(run_tree, 'LHEScaleSumw[%d]*genEventSumw' % i)
                                        for i in range(nScaleWeights)])
        print('LHEScaleWeightNorm: ' + str(scale_weight_norm_buff))
        _fill_const_branch(tree, "LHEScaleWeightNorm", scale_weight_norm_buff, lenVar=nScaleWeights)
        
    if tree.GetBranch('LHEPdfWeight'):
        run_tree.GetEntry(0)
        nPdfWeights = run_tree.nLHEPdfSumw
        pdf_weight_norm_buff = array('f',
                                     [sumwgts / _get_sum(run_tree, 'LHEPdfSumw[%d]*genEventSumw' % i)
                                      for i in range(nPdfWeights)])
        print('LHEPdfWeightNorm: ' + str(pdf_weight_norm_buff))
        _fill_const_branch(tree, "LHEPdfWeightNorm", pdf_weight_norm_buff, lenVar=nScaleWeights)

    # fill PS weight re-normalization factors
    if tree.GetBranch('PSWeight') and run_tree.GetBranch('PSSumw'):
        run_tree.GetEntry(0)
        nPSWeights = run_tree.nPSSumw
        ps_weight_norm_buff = array('f',
                                    [sumwgts / _get_sum(run_tree, 'PSSumw[%d]*genEventSumw' % i)
                                     for i in range(nPSWeights)])
        print('PSWeightNorm: ' + str(ps_weight_norm_buff))
        _fill_const_branch(tree, 'PSWeightNorm', ps_weight_norm_buff, lenVar='nPSWeight')
               
    tree.Write(treename, ROOT.TObject.kOverwrite)
    f.Close()

def run_add_weights():
    xsec_dict = parse_sample_xsec(args.xsec_file)

    with open(args.jobs_dir + "/metadata.json", 'r') as file:
        data = json.load(file) 
    
    year = data["year"]
    dataset_type = data["type"]
    base_output_dir = data["output_dir"]
    
    sample_dirs = [d.name for d in Path(os.path.join(base_output_dir, dataset_type, year)).iterdir() if d.is_dir()]

    for sample in sample_dirs:
        physics_process_dirs = [d.name for d in Path(os.path.join(base_output_dir, dataset_type, year, sample)).iterdir() if d.is_dir()]

        for physics_process in physics_process_dirs:
            if physics_process not in xsec_dict: 
                print(f"Process {physics_process} not found in xsec file, xsec not added") 
                continue
                     
            xsec = xsec_dict[physics_process]
            print(f"Adding weights for physics process {physics_process}: xsec = {xsec}")  
            
            for file in Path(os.path.join(base_output_dir, dataset_type, year, sample, physics_process)).iterdir():
                add_weights(file, xsec)

def check_job_status():
    
    file_path = args.jobs_dir + '/metadata.json'
    with open(file_path, 'r') as file:
        data = json.load(file)
    
    njobs = len(data['jobs'])
    jobids = {'running': [], 'failed': [], 'completed': []}
    for jobid in range(njobs):
        logpath = os.path.join(args.jobs_dir, "log", '%d.log' % jobid)
        # print(logpath)
        if not os.path.exists(logpath):
            print('Cannot find log file %s' % logpath)
            jobids['failed'].append(str(jobid))
            continue
        with open(logpath) as logfile:
            errormsg = None
            finished = False
            for line in reversed(logfile.readlines()):
                if 'Job removed' in line or 'aborted' in line:
                    errormsg = line
                if 'Job submitted from host' in line:
                    # if seeing this first: the job has been resubmited
                    break
                if 'return value' in line:
                    if 'return value 0' in line:
                        finished = True
                    else:
                        errormsg = line
                    break
            if errormsg:
                print(logpath + '\n   ' + errormsg)
                jobids['failed'].append(str(jobid))
            else:
                if finished:
                    jobids['completed'].append(str(jobid))
                else:
                    jobids['running'].append(str(jobid))
    assert sum(len(jobids[k]) for k in jobids) == njobs
    all_completed = len(jobids['completed']) == njobs
    info = {k: len(jobids[k]) for k in jobids if len(jobids[k])}
    print('Job %s status: ' % args.jobs_dir + str(info))
    return all_completed, jobids

def resubmit():
    
    jobids = check_job_status(args)[1]['failed']
    jobids_file = os.path.join(args.jobs_dir, 'resubmit.txt')

    with open(jobids_file, 'w') as f:
        f.write('\n'.join(jobids))
    
    write_condor_submit(jobids_file="resubmit.txt")
        
def main():

    jobs_dir = "jobs_" + args.type + "_" + args.year
    args.jobs_dir = jobs_dir

    if args.resubmit:
        resubmit()
        sys.exit(0)
    
    if args.check_status:
        check_job_status()
        sys.exit(0)
   
    if args.post:
        # run_add_weights()
        merge_output_files()
        sys.exit(0)   
        
    helpers.check_if_dir_exists(jobs_dir)
    helpers.check_if_dir_exists(args.output)

    print("Will write output trees to " + args.output)
    print("Number of files per job: " + str(args.n))
    
    ## create necessary dirs
    os.system("mkdir -p " + jobs_dir + "/log/")
    
    ## copy necessary files to jobs dir
    os.system("cp static_files/processor.py " + jobs_dir)
    os.system("cp static_files/keep_and_drop*.txt " + jobs_dir)
    os.system("cp static_files/condor_exec.sh " + jobs_dir)
    
    ## create metadata json and condor submit files in the jobs dir
    create_metadata_json()
    create_condor_submit()

if __name__ == "__main__":
    main()