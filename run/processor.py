from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoHc.producers.HcTreeProducer import HcTreeProducerModule

import sys
import json 

jobid = name = sys.argv[1]

## read json file
file_path = 'metadata.json'
with open(file_path, 'r') as file:
    data = json.load(file)
    
files = data["jobs"][int(jobid)]["input_files"]
output_dir = data["output_dir"]
sample = data["jobs"][int(jobid)]["sample_name"]

## convert keep_and_drop.txt to python list
with open('keep_and_drop.txt', 'r') as file:
    keep_and_drop_branches = file.readlines()
keep_and_drop_branches = [line.strip() for line in keep_and_drop_branches]


p = PostProcessor(
    outputDir = output_dir, 
    inputFiles = files, 
    modules=[HcTreeProducerModule()],
    outputbranchsel=keep_and_drop_branches,
    postfix="_"+sample,
    prefetch=True)
p.run()