import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}

class JetJERC(Module, object):
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

            self.out.branch('jetJERC', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        # FIX ME 
        
        return True
    
class JetVMAP(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jetvetomaps.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)

    def get_sf(self, sf_type, lep):

        return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('jetVMAP', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        # FIX ME 
        
        return True