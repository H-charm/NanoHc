import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}
key_dict = {"2022":     'Summer22_23Sep2023_RunCD_V1',
            "2022EE":   'Summer22EE_23Sep2023_RunEFG_V1',
            "2023":     'Summer23Prompt23_RunC_V1',
            "2023BPix": 'Summer23BPixPrompt23_RunD_V1'
            }

class JetJERCProducer(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]

        correction_file_JERC = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jet_jerc.json.gz'
        correction_file_JERsmear = '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz'

        self.corr_JERC = correctionlib.CorrectionSet.from_file(correction_file_JERC)
        self.corr_JERsmear = correctionlib.CorrectionSet.from_file(correction_file_JERsmear)
    
    def get_sf(self, sf_type, lep):
        
        return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('jetJERCWeight', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        # FIX ME 
        
        return True
    
class JetVMAPProducer(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jetvetomaps.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)[key_dict[self.year]]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('jetVMAPWeight', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # According to JetMet POG: These maps should be applied similarly both on Data and MC, to keep the phase-spaces equal.
        # https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps

        for jet in event.selectedJets:

            wgtjetVMAP = self.corr.evaluate("jetvetomap", jet.eta, jet.phi)

        self.out.fillBranch('jetVMAPWeight', wgtjetVMAP)
        
        return True