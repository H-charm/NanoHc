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

"""
Add jetId variable, following the example in https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py .
See example in test/example_jetId.py for usage.
"""
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
from array import array

class JetIdProducer(Module):
    def __init__(self, year, dataset_type,  jetType="AK4PUPPI"):
        """Module to determine jetID variables (passTight, passTightLepVeto), 
        packed in Jet_jetId as in nanoAODv12
        Parameters:
        - json: jetID json file
        - jetType: "AK4PUPPI" or "AK4CHS"

        #TODO THIS SELECTION WAS ON POSTPROCESSOR MAYBE INCLUDE IT nJet>0 && Jet_pt>30? 
        """

        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jetid.json.gz'
        self.evaluator = correctionlib.CorrectionSet.from_file(correction_file)
        self.key_tight = f"{jetType}_Tight"
        self.key_tightLeptonVeto = f"{jetType}_TightLeptonVeto"
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_jetId", "b", lenVar="nJet", title="Jet ID flag: bit2 is tight, bit3 is tightLepVeto (recomputed using JSON)") # Save as UChar_t as it was up to nanoAODv14


    def analyze(self, event):
        jets = Collection(event, "Jet")

        jet_Ids = array('B', event.nJet*[0]) # Note: UChar_t is uppercase 'B' in python array
        for ijet, jet in enumerate(jets):
            multiplicity = jet.chMultiplicity + jet.neMultiplicity

            passTight = self.evaluator[self.key_tight].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )

            passTightLepVeto = self.evaluator[self.key_tightLeptonVeto].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )
            
            jet_Ids[ijet] = int(passTight)*2 + int(passTightLepVeto)*4

        self.out.fillBranch("Jet_jetId", jet_Ids)
        return True