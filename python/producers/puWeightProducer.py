import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}
key_dict = {"2022":     'Collisions2022_355100_357900_eraBCD_GoldenJson',
            "2022EE":   'Collisions2022_359022_362760_eraEFG_GoldenJson',
            "2023":     'Collisions2023_366403_369802_eraBC_GoldenJson',
            "2023BPix": 'Collisions2023_369803_370790_eraD_GoldenJson'
            }


class PileupWeightProducer(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.nvtxVar = 'Pileup_nTrueInt'
        self.era = era_dict[self.year]
        correction_file = f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{self.era}/puWeights.json.gz"
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)[key_dict[self.year]]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('puWeight', "F")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if not self.isMC:
            return True
        
        puWeight = self.corr.evaluate(getattr(event, self.nvtxVar), "nominal")

        self.out.fillBranch('puWeight', puWeight)

        return True