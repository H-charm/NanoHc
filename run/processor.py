from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoHc.producers.HcTreeProducer import HcTreeProducer
from PhysicsTools.NanoHc.producers.leptonvariables import LeptonVariablesModule
from PhysicsTools.NanoHc.producers.topleptonmva import TopLeptonMvaModule

import sys
import json 

jobid = name = sys.argv[1]

## read json file
file_path = 'metadata.json'
with open(file_path, 'r') as file:
    data = json.load(file)
    
files = data["jobs"][int(jobid)]["input_files"]
output_dir = data["output_dir"]
year = data["year"]
sample = data["jobs"][int(jobid)]["sample_name"]

## convert keep_and_drop_input.txt to python list
with open('keep_and_drop_input.txt', 'r') as file:
    keep_and_drop_input_branches = file.readlines()
keep_and_drop_input_branches = [line.strip() for line in keep_and_drop_input_branches]

## convert keep_and_drop_output.txt to python list
with open('keep_and_drop_output.txt', 'r') as file:
    keep_and_drop_output_branches = file.readlines()
keep_and_drop_output_branches = [line.strip() for line in keep_and_drop_output_branches]

p = PostProcessor(
    outputDir = output_dir, 
    inputFiles = files, 
    modules=[HcTreeProducer(year),
             LeptonVariablesModule(),
             TopLeptonMvaModule("2018", 'ULv1')],
    branchsel=keep_and_drop_input_branches,
    outputbranchsel=keep_and_drop_output_branches,
    postfix="_"+sample,
    prefetch=True)
p.run()