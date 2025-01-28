from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoHc.producers.BaselineProducer import BaselineProducer
from PhysicsTools.NanoHc.producers.puWeightProducer import PileupWeightProducer
from PhysicsTools.NanoHc.producers.leptonSFProducer import ElectronSFProducer, MuonSFProducer
from PhysicsTools.NanoHc.producers.leptonvariables import LeptonVariablesModule
from PhysicsTools.NanoHc.producers.topleptonmva import TopLeptonMvaModule
from PhysicsTools.NanoHc.producers.jetSFProducer import JetVMAPProducer

import sys
import json 
import os

jobid = name = sys.argv[1]

## read json file
file_path = 'metadata.json'
with open(file_path, 'r') as file:
    data = json.load(file)

## load from json file    
files = data["jobs"][int(jobid)]["input_files"]
base_output_dir = data["output_dir"]
year = data["year"]
dataset_type = data["type"]
golden_json = data["golden_json"]
sample = data["jobs"][int(jobid)]["sample_name"]
physics_process = data["jobs"][int(jobid)]["physics_process"]

## convert keep_and_drop_input.txt to python list
with open('keep_and_drop_input.txt', 'r') as file:
    keep_and_drop_input_branches = file.readlines()
keep_and_drop_input_branches = [line.strip() for line in keep_and_drop_input_branches]

## convert keep_and_drop_output.txt to python list
with open('keep_and_drop_output.txt', 'r') as file:
    keep_and_drop_output_branches = file.readlines()
keep_and_drop_output_branches = [line.strip() for line in keep_and_drop_output_branches]

output_dir = os.path.join(base_output_dir, dataset_type, year, sample, physics_process) 

p = PostProcessor(
    outputDir = output_dir, 
    inputFiles = files, 
    modules=[
            # LeptonVariablesModule(),
            # TopLeptonMvaModule(year, 'ULv2'),     
            BaselineProducer(year, dataset_type, sample),
            JetVMAPProducer(year,dataset_type),
            # PileupWeightProducer(year, dataset_type),
            # ElectronSFProducer(year, dataset_type), # pt binning starts at 10, our selections at 7 (keep it out for now)
            # MuonSFProducer(year, dataset_type),
            ],
    branchsel=keep_and_drop_input_branches,
    outputbranchsel=keep_and_drop_output_branches,
    postfix="_" + physics_process,
    prefetch=True,
    jsonInput=golden_json
    )
p.run()