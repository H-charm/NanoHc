from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoHc.producers.BaselineProducer import BaselineProducer
from PhysicsTools.NanoHc.producers.jetSFProducer import JetVMAPProducer 
from PhysicsTools.NanoHc.producers.jetIDProducer import jetId 
# from PhysicsTools.NanoHc.producers.puWeightProducer import PileupWeightProducer
# from PhysicsTools.NanoHc.producers.leptonSFProducer import ElectronSFProducer, MuonSFProducer
# from PhysicsTools.NanoHc.producers.jetSFProducer import JetVMAPProducer, jetJERCProducer
# from PhysicsTools.NanoHc.producers.leptonScaleResProducer import eleScaleRes, muonScaleRes

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
if dataset_type == "data":
    era_data = data["jobs"][int(jobid)]["era"]
else:
    era_data = None

# Choose file names based on dataset type
if dataset_type == "data":
    input_file = "keep_and_drop_input_data.txt"
    output_file = "keep_and_drop_output_data.txt"
else:
    input_file = "keep_and_drop_input_mc.txt"
    output_file = "keep_and_drop_output_mc.txt"

# Convert input file to python list
with open(input_file, 'r') as file:
    keep_and_drop_input_branches = file.readlines()
keep_and_drop_input_branches = [line.strip() for line in keep_and_drop_input_branches]

# Convert output_file to python list
with open(output_file, 'r') as file:
    keep_and_drop_output_branches = file.readlines()
keep_and_drop_output_branches = [line.strip() for line in keep_and_drop_output_branches]

output_dir = os.path.join(base_output_dir, dataset_type, year, sample, physics_process) 

p = PostProcessor(
    outputDir = output_dir, 
    inputFiles = files, 
    modules=[
            JetVMAPProducer(year,dataset_type),
            jetId(year,dataset_type),
            # jetJERCProducer(year, era_data, dataset_type),
            # eleScaleRes(year,dataset_type),
            # muonScaleRes(year,dataset_type),
            BaselineProducer(year, dataset_type, sample),
            # PileupWeightProducer(year, dataset_type),
            # ElectronSFProducer(year, dataset_type), # pt binning starts at 10, our selections at 7 -- For 2022 manually fixed, for 2023 it passes as 1. 
            # MuonSFProducer(year, dataset_type),
            ],
    branchsel=keep_and_drop_input_branches,
    outputbranchsel=keep_and_drop_output_branches,
    postfix="_" + physics_process,
    prefetch=True,
    jsonInput=golden_json
    )
p.run()