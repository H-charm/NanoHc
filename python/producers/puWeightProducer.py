import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


era_dict = {"2021": '2022_Summer22', "2022": '2022_Summer22EE'}
key_dict = {"2021": 'Cert_Collisions2022_355100_362760_Golden.json',
            "2022": 'Cert_Collisions2022_355100_362760_Golden.json'
            }


class PileupWeightProducer(Module, object):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.nvtxVar = 'Pileup_nTrueInt'
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{self.era}/puWeights.json.gz'
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