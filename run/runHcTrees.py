#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import yaml
import json
import os
import sys
from pathlib import Path

golden_json = {
    '2015': 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    '2016': 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    '2017': 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
    '2018': 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',
}

## to divide datasets by number of files per job
def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


## read samples yaml file and produce json file to be used by condor
def create_metadata_json(args):
    
    dataset_type = args.type
    jobs_dir_name = "jobs_" + dataset_type + "_" + args.year
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
    json_file = jobs_dir_name + '/metadata.json'
    json_content = {}
    
    json_content["output_dir"] = args.output
    json_content["year"] = args.year    
    json_content["type"] = dataset_type
    if dataset_type == "data":  
        json_content["golden_json"] = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoHc/data/JSON/" + golden_json[args.year]
    else:
        json_content["golden_json"] = None
    json_content["sample_names"] = []
    json_content["physics_processes"] = []
    json_content["jobs"] = []
    job_id = 0
    # json_content["sample_names"].append([sample_name for sample_name in samples])
    
    for sample_name in samples: json_content["sample_names"].append(sample_name)
    for physics_process in physics_processes: json_content["physics_processes"].append(physics_process)
    
    for sample in samples:
        for physics_process in das_dict[sample]:
            for chunk in enumerate(get_chunks(das_dict[sample][physics_process],args.n)):
                json_content["jobs"].append({"job_id": job_id, "input_files": chunk[1], "sample_name": sample ,"physics_process": physics_process})
                job_id += 1

    with open(json_file, 'w') as file:
        json.dump(json_content, file, indent=4)

## produce condor submit file
def create_condor_submit(jobs_dir_name):
    
    cmssw_base = os.environ['CMSSW_BASE']
    jobs_dir_path = os.getcwd() + "/" + jobs_dir_name

    ## find number of jobs from json file
    with open(jobs_dir_name + "/metadata.json", 'r') as file:
        data = json.load(file)  
    njobs = len(data["jobs"])
    jobs_list = [i for i in range(njobs)]
    
    ## write jobs txt file
    with open(jobs_dir_name + "/job_ids.txt", 'w') as file:
        for index, item in enumerate(jobs_list):
            if index < len(jobs_list) - 1:
                file.write(str(item) + '\n')
            else:
                file.write(str(item))
        
    ## write condor submit 
    condor_submit_file = open(jobs_dir_name + "/submit.sh","w")
    condor_submit_file.write('''
executable = condor_exec.sh

arguments = $(jobid) ''' + cmssw_base + ''' ''' + jobs_dir_path + ''' 

request_memory  = 2000
request_disk    = 10000000

error   = log/err.$(Process)
output  = log/out.$(Process)
log     = log/logFile.log

+JobFlavour = "longlunch"

queue jobid from job_ids.txt
    ''')

    condor_submit_file.close()


def merge_output_files(jobs_dir_name):

    file_path = jobs_dir_name + '/metadata.json'
    with open(file_path, 'r') as file:
        data = json.load(file)
        
    base_output_dir = data["output_dir"]
    dataset_type = data["type"]
    year = data["year"]
    sample_names = data["sample_names"]
    merged_dir = os.path.join(base_output_dir, dataset_type, year, "merged")
    os.system("mkdir -p " + merged_dir)
    
    for sample_name in sample_names:
        input_files = []
        output_file = os.path.join(merged_dir, sample_name + "_tree.root")
        
        physics_process_dirs = [d.name for d in Path(os.path.join(base_output_dir, dataset_type, year, sample_name)).iterdir() if d.is_dir()]
        for physics_process_dir in physics_process_dirs:
            infiles_dir = os.path.join(base_output_dir, dataset_type, year, sample_name, physics_process_dir)
            
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

def add_weight_branch(file, xsec, lumi=1000., treename='Events', wgtbranch='xsecWeight'):
    from array import array
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    def _get_sum(tree, wgtvar):
        htmp = ROOT.TH1D('htmp', 'htmp', 1, 0, 10)
        tree.Project('htmp', '1.0', wgtvar)
        return float(htmp.Integral())

    def _fill_const_branch(tree, branch_name, buff, lenVar=None):
        if lenVar is not None:
            b = tree.Branch(branch_name, buff, '%s[%s]/F' % (branch_name, lenVar))
            b_lenVar = tree.GetBranch(lenVar)
            buff_lenVar = array('I', [0])
            b_lenVar.SetAddress(buff_lenVar)
        else:
            b = tree.Branch(branch_name, buff, branch_name + '/F')

        b.SetBasketSize(tree.GetEntries() * 2)  # be sure we do not trigger flushing
        for i in range(tree.GetEntries()):
            if lenVar is not None:
                b_lenVar.GetEntry(i)
            b.Fill()

        b.ResetAddress()
        if lenVar is not None:
            b_lenVar.ResetAddress()
    
    f = ROOT.TFile(str(file), 'UPDATE')
    run_tree = f.Get('Runs')
    tree = f.Get(treename)

    # fill cross section weights to the 'Events' tree
    sumwgts = _get_sum(run_tree, 'genEventSumw')
    xsecwgt = xsec * lumi / sumwgts
    xsec_buff = array('f', [xsecwgt])
    _fill_const_branch(tree, wgtbranch, xsec_buff)

    tree.Write(treename, ROOT.TObject.kOverwrite)
    f.Close()


def run_lumi_wgt(args, jobs_dir_name):
    xsec_dict = parse_sample_xsec(args.xsec_file)

    with open(jobs_dir_name + "/metadata.json", 'r') as file:
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
            print(f"Adding lumi wgt for physics process {physics_process}: xsec = {xsec}")  
            
            for file in Path(os.path.join(base_output_dir, dataset_type, year, sample, physics_process)).iterdir():
                add_weight_branch(file, xsec)


def main():

    ## parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=str, help='Year to run', default="2018")
    parser.add_argument('--output', type=str, help='Output dir', default = "/eos/user/i/iparaske/HcTrees/")
    parser.add_argument('--type', type=str, help='mc or data', default = "mc", choices=['mc', 'data'])
    parser.add_argument('--post',help='Merge output files',action='store_true')
    parser.add_argument('-n',type=int, help='Number of files per job', default=10)
    parser.add_argument('--xsec-file', type=str, help='xsec file', default = "samples/xsec.conf")
    args = parser.parse_args()
    
    jobs_dir_name = "jobs_" + args.type + "_" + args.year

    if args.post:
        run_lumi_wgt(args, jobs_dir_name)
        merge_output_files(jobs_dir_name)
        sys.exit(0)

    ## print info    
    print("Will write output trees to " + args.output)
    print("Number of files per job: " + str(args.n))
    
    ## create jobs dir based on year
    os.system("mkdir -p " + jobs_dir_name)

    ## create logs dir
    os.system("mkdir -p " + jobs_dir_name + "/log/")
    
    ## copy processor to jobs dir
    os.system("cp processor.py " + jobs_dir_name)
    
    ## copy branch selection to jobs dir
    os.system("cp keep_and_drop_input.txt " + jobs_dir_name)
    os.system("cp keep_and_drop_output.txt " + jobs_dir_name)
    
    ## copy condor executable to jobs dir
    os.system("cp condor_exec.sh " + jobs_dir_name)
    
    ## create metadata json and condor submit files in the jobs dir
    create_metadata_json(args)
    create_condor_submit(jobs_dir_name)

if __name__ == "__main__":
    main()