import os

# Set up the NATModules muonScaleRes module
def getMuonScaleRes(year, dataset_type, overwritePt=True) :
    from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes 

    key_dict={
        "2022":     "2022_Summer22.json.gz",
        "2022EE":   "2022_Summer22EE.json.gz",
        "2023":     "2023_Summer23.json.gz",
        "2023BPix": "2023_Summer23BPix.json.gz"
    }

    is_mc = True if dataset_type == "mc" else False

    # Json files for Muons Scale and Smearing corrections are taken from https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/tree/master/corrections
    json = f"../../data/MuonScale/{key_dict[year]}" 

    print("***muonScaleRes: year:", year, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", json)
    return muonScaleRes(json, is_mc, overwritePt, minPt=3.)
