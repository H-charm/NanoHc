import json
import gzip

def fix_json_gz(input_gz_file, output_gz_file):
    # Step 1: Read & Decompress the .json.gz file
    with gzip.open(input_gz_file, "rt", encoding="utf-8") as f:
        data = json.load(f)

    # Step 2: Recursive function to replace "inf" and "-inf"
    def replace_inf(obj):
        if isinstance(obj, dict):
            return {k: replace_inf(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_inf(v) for v in obj]
        elif obj == "inf":
            return float("Infinity")
        elif obj == "-inf":
            return float("-Infinity")
        return obj

    # Step 3: Add '1' to the "run" bin edges
    # def add_run_edge(obj):
    #     if isinstance(obj, dict):
    #         for key, value in obj.items():
    #             obj[key] = add_run_edge(value)
    #     elif isinstance(obj, list):
    #         for i in range(len(obj)):
    #             obj[i] = add_run_edge(obj[i])
    #     elif isinstance(obj, float) or isinstance(obj, int):  # Ensure numerical check
    #         return obj
        
    #     # Check if "edges" exist inside a "nodetype" multibinning with "run"
    #     if isinstance(obj, dict) and obj.get("nodetype") == "multibinning" and "run" in obj.get("inputs", []):
    #         for edge_list in obj.get("edges", []):
    #             if isinstance(edge_list, list) and 1 not in edge_list:
    #                 edge_list.append(1)  # Add '1' if missing
    #                 edge_list.sort()  # Ensure it remains sorted
    #     return obj

    # Apply both modifications
    fixed_data = replace_inf(data)
    # fixed_data = add_run_edge(fixed_data)

    # Step 4: Recompress & Save the fixed JSON as a .json.gz file
    with gzip.open(output_gz_file, "wt", encoding="utf-8") as f:
        json.dump(fixed_data, f, indent=4)

    print(f"✅ Fixed .json.gz file saved as {output_gz_file}")

# Usage Example
input_gz_json = "/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/ElectronScale/electronSS_preBPix.json.gz"
output_gz_json = "/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/ElectronScale/electronSS_preBPix_inf.json.gz"


fix_json_gz(input_gz_json, output_gz_json)
