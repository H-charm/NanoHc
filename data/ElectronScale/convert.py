import json
import gzip

def fix_multibinning(input_gz_file, output_gz_file):
    # Step 1: Read & decompress the JSON file
    with gzip.open(input_gz_file, "rt", encoding="utf-8") as f:
        data = json.load(f)

    # Step 2: Navigate to the "corrections" section
    for correction in data.get("corrections", []):
        if correction.get("name") == "2023PromptC_ScaleJSON":
            # Navigate through nested categories
            for content in correction["data"]["content"]:
                if content["value"]["nodetype"] == "category":
                    for sub_content in content["value"]["content"]:
                        if sub_content["value"]["nodetype"] == "multibinning":
                            # Fix MultiBinning edges & content
                            for i, edge_list in enumerate(sub_content["value"]["edges"]):
                                content_length = len(sub_content["value"]["content"])
                                edge_length = len(edge_list)

                                # Ensure bin edges are one more than content nodes
                                if edge_length != content_length + 1:
                                    print(f"⚠️ Fixing MultiBinning mismatch in {sub_content['value']['inputs']}")
                                    if edge_length > content_length + 1:
                                        sub_content["value"]["edges"][i] = edge_list[:content_length + 1]  # Trim extra edges
                                    else:
                                        while len(sub_content["value"]["edges"][i]) < content_length + 1:
                                            sub_content["value"]["edges"][i].append(sub_content["value"]["edges"][i][-1] + 1)  # Add missing edges
                                
                                # Ensure "run" edges include 1 and remain sorted
                                if "run" in sub_content["value"]["inputs"]:
                                    if 1 not in sub_content["value"]["edges"][i]:
                                        sub_content["value"]["edges"][i].append(1)
                                        sub_content["value"]["edges"][i].sort()

    # Step 3: Save the fixed JSON file
    with gzip.open(output_gz_file, "wt", encoding="utf-8") as f:
        json.dump(data, f, indent=4)

    print(f"✅ Fixed .json.gz file saved as {output_gz_file}")

# Usage Example
input_gz_json = "/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/ElectronScale/electronSS_preBPix_mod.json.gz"
output_gz_json = "/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/ElectronScale/electronSS_preBPix_modi.json.gz"

fix_multibinning(input_gz_json, output_gz_json)
