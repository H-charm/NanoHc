import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}
key_dict = {
    "2022":     'Summer22_23Sep2023_RunCD_V1',
    "2022EE":   'Summer22EE_23Sep2023_RunEFG_V1',
    "2023":     'Summer23Prompt23_RunC_V1',
    "2023BPix": 'Summer23BPixPrompt23_RunD_V1'
}
    
class JetVMAPProducer(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        """According to JetMet POG:
           These are the jet veto maps showing regions with an excess of jets (hot zones) and lack of jets
           (cold zones). Using the phi-symmetry of the CMS detector, these areas with detector and or calibration issues can be pinpointed.
           This module applies veto maps that veto out events with important jets in these "hot" or "cold" zones. 
           These maps should be applied similarly both on Data and MC, to keep the phase-spaces equal.
           Non-zero value indicates that the region is vetoed 
           https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
        """

        self.year = year
        self.dataset_type = dataset_type
        self.veto_map_name = "jetvetomap"
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jetvetomaps.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)[key_dict[self.year]]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Flag_JetVetoed", "O", title="Event veto flag from Jet Veto Map")

    def fixPhi(self, phi):
        epsilon = 1e-6  # Small offset to avoid boundary issues
        if phi > np.pi:
            #print(f"phi {phi} is greater than pi. Setting phi to pi - epsilon.")
            phi = np.pi - epsilon
        elif phi < -np.pi:
            #print(f"phi {phi} is less than -pi. Setting phi to -pi + epsilon.")
            phi = -np.pi + epsilon
        return phi

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)

        nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        - jets that don't overlap with PF muon (dR < 0.2)
        """
        jets = Collection(event, "Jet")
        veto_flag = False 
        for i,jet in enumerate(jets):
            if (jet.pt> 15 and (jet.jetId ==2 or jet.jetId ==6) and jet.chEmEF < 0.9 and jet.neEmEF <0.9 and jet.muonIdx1 == -1 and jet.muonIdx2 == -1):
            
                # Correct phi and evaluate veto map
                phi = self.fixPhi(jet.phi)
                veto_map_value = self.corr.evaluate(self.veto_map_name, jet.eta, phi)

                # Check if the jet is vetoed
                if veto_map_value > 0:
                    veto_flag = True  # Set flag if a vetoed jet is found
                    break  # Break out of the loop since we only need one veto to trigger

        # Fill the branch with the veto result
        self.out.fillBranch("Flag_JetVetoed", veto_flag)
        
        return True