import numpy as np
import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

era_dict = {
    "2022": '2022',
    "2022EE": '2022_EE',
    "2023": '2023',
    "2023BPix": '2023_BPix'
}

class MuonTriggerProducer(Module):
    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]

        self.root_file = ROOT.TFile(f'../../data/MuonTRGSF/{self.year}/ScaleFactors_Muon_Z_HLT_{self.era}_eta_pt.root', "READ")
        #self.h_Mu_SF = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt").Clone()   
        self.h_Mu_MCEff = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyMC").Clone()
        self.h_Mu_DataEff = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyData").Clone()
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = self.dataset_type == "mc"
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('muTriggerWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muTriggerWeightUp', "F", limitedPrecision=10)
                self.out.branch('muTriggerWeightDown', "F", limitedPrecision=10)
    
    def getSF(self, pt, eta, corr_type):
        if corr_type == 'mc':
            h = self.h_Mu_MCEff
        else:
            h = self.h_Mu_DataEff
        # h = self.h_Mu_SF
        binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
        binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
        return h.GetBinContent(binX, binY)

    def analyze(self, event):
        if not self.isMC:
            return True

        sf_list = []
            
        # if len(event.selectedMuons) == 2:
        #     mu1, mu2 = event.selectedMuons[:2]
        if len(event.fullIDMuons) == 2:
            mu1, mu2 = event.fullIDMuons[:2]

            if abs(mu1.pdgId) != 13 or abs(mu2.pdgId) !=13:
                return True
            
            pt1 = min(mu1.pt, 199.0)
            pt2 = min(mu2.pt, 199.0)

            eta1 = mu1.eta
            eta2 = mu2.eta

            if pt1 < 26.0:
                effMC1 = 1.0
                effData1 = 1.0
            else:
                effMC1 = self.getSF(pt1, eta1, 'mc')
                effData1 = self.getSF(pt1, eta1, 'data')
            
            if pt2 < 26.0:
                effMC2 = 1.0
                effData2 =1.0
            else:
                effMC2 = self.getSF(pt2, eta2, 'mc')
                effData2 = self.getSF(pt2, eta2, 'data')

            SF_den = effMC1 + effMC2 - effMC1 * effMC2

            SF_num = effData1 + effData2 - effData1 * effData2

            SF = SF_num / SF_den

            sf_list.append(SF)

        HLT_weight = np.prod(sf_list) if sf_list else 1.0

        self.out.fillBranch("muTriggerWeight", HLT_weight)

        return True
