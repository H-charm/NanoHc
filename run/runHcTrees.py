#!/usr/bin/env python

from os.path import samefile, sameopenfile
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import yaml
import json
import os
import sys
import subprocess

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

    das_files = {}
    for sample in samples:
        for sub_sample in samples[sample]:
            ## find files using DAS
            print("DAS query for dataset " + sub_sample)
            das_query = 'dasgoclient --query="file dataset=' + sub_sample + '"'
            query_out = os.popen(das_query)
            files_found = ['root://xrootd-cms.infn.it/'+_file.strip() for _file in query_out]
            das_files.setdefault(sample, []).extend(files_found)
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
    json_content["jobs"] = []
    job_id = 0
    for sample_name in samples:
        json_content["sample_names"].append(sample_name)
        for chunk in enumerate(get_chunks(das_files[sample_name],args.n)):
            json_content["jobs"].append({"job_id": job_id, "input_files": chunk[1], "sample_name": sample_name})
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


def merge_output_files(output_dir, jobs_dir_name):

    file_path = jobs_dir_name + '/metadata.json'
    with open(file_path, 'r') as file:
        data = json.load(file)
    output_dir = data["output_dir"]
    sample_names = data["sample_names"]

    for sample_name in sample_names:
        input_files = []
        output_file = output_dir + sample_name + "_tree.root"
        for root, _, files in os.walk(output_dir):
            for file in files:
                if '_' in file:
                    _, second_part = file.split('_', 1)
                    if second_part.startswith(sample_name):
                        input_files.append(os.path.join(root, file))
                        
        cmd = 'haddnano.py {outfile} {infiles}'.format(outfile=output_file, infiles=' '.join(input_files))
        os.system(cmd)

def main():

    ## parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=str, help='Year to run', default="2018")
    parser.add_argument('--output', type=str, help='Output dir', default = "/eos/user/i/iparaske/HcTrees/")
    parser.add_argument('--type', type=str, help='mc or data', default = "mc", choices=['mc', 'data'])
    parser.add_argument('--post',help='Merge output files',action='store_true')
    parser.add_argument('-n',type=int, help='Number of files per job', default=10)
    args = parser.parse_args()
    
    jobs_dir_name = "jobs_" + args.type + "_" + args.year

    if args.post:
        merge_output_files(args.output, jobs_dir_name)
        sys.exit(0)

    ## print info    
    print("Will write output trees to " + args.output)
    print("Number of files per job: " + str(args.n))
    
    ## create output dir
    os.system("mkdir -p " + args.output)
    
    ## create jobs dir based on year
    os.system("mkdir -p " + jobs_dir_name)

    ## create logs dir
    os.system("mkdir -p " + jobs_dir_name + "/log/")
    
    ## copy processor to jobs dir
    os.system("cp processor.py " + jobs_dir_name)
    
    ## copy branches selections to jobs dir
    os.system("cp keep_and_drop_input.txt " + jobs_dir_name)
    os.system("cp keep_and_drop_output.txt " + jobs_dir_name)
    
    ## copy condor executable to jobs dir
    os.system("cp condor_exec.sh " + jobs_dir_name)
    
    ## create metadata json and condor submit files in the jobs dir
    create_metadata_json(args)
    create_condor_submit(jobs_dir_name)

if __name__ == "__main__":
    main()