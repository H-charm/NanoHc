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
        self.h_Mu_MCEff_Up = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyMC_errU").Clone()
        self.h_Mu_MCEff_Down = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyMC_errD").Clone()
        self.h_Mu_DataEff = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyData").Clone()
        self.h_Mu_DataEff_Up = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyData_errU").Clone()
        self.h_Mu_DataEff_Down = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyData_errD").Clone()
        
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
            if self.doSysVar:
                h_up = self.h_Mu_MCEff_Up
                h_down = self.h_Mu_MCEff_Down
        else:
            h = self.h_Mu_DataEff
            if self.doSysVar:
                h_up = self.h_Mu_DataEff_Up
                h_down = self.h_Mu_DataEff_Down
        # h = self.h_Mu_SF
        if self.doSysVar:
            binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
            binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())

            binX_Up = min(max(1, h_up.GetXaxis().FindBin(eta)), h_up.GetNbinsX())
            binY_Up = min(max(1, h_up.GetYaxis().FindBin(pt)), h_up.GetNbinsY())
            binX_Down = min(max(1, h_down.GetXaxis().FindBin(eta)), h_down.GetNbinsX())
            binY_Down = min(max(1, h_down.GetYaxis().FindBin(pt)), h_down.GetNbinsY())
            
            return h.GetBinContent(binX, binY), h_up.GetBinContent(binX_Up, binY_Up), h_down.GetBinContent(binX_Down, binY_Down)
        else:
            binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
            binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())

            return h.GetBinContent(binX, binY)

    def analyze(self, event):
        if not self.isMC:
            return True

        if self.doSysVar:
            sf_list = []
            sf_list_up = []
            sf_list_down = []

            if len(event.selectedMuons) == 2:
                mu1, mu2 = event.selectedMuons[:2]

                if abs(mu1.pdgId) != 13 or abs(mu2.pdgId) !=13:
                    return True
                
                pt1 = min(mu1.pt, 199.0)
                pt2 = min(mu2.pt, 199.0)

                eta1 = mu1.eta
                eta2 = mu2.eta

                if pt1 < 26.0:
                    effMC1 = 1.0 
                    effMC1_up = 1.0 
                    effMC1_down = 1.0
                    effData1 = 1.0
                    effData1_up = 1.0
                    effData1_down = 1.0
                else:
                    effMC1, effMC1_up, effMC1_down = self.getSF(pt1, eta1, 'mc')
                    effData1, effData1_up, effData1_down = self.getSF(pt1, eta1, 'data')
                
                if pt2 < 26.0:
                    effMC2 = 1.0
                    effMC2_up = 1.0
                    effMC2_down = 1.0
                    effData2 =1.0
                    effData2_up = 1.0
                    effData2_down = 1.0
                else:
                    effMC2, effMC2_up, effMC2_down = self.getSF(pt2, eta2, 'mc')
                    effData2, effData2_up, effData2_down = self.getSF(pt2, eta2, 'data')

                SF_den = effMC1 + effMC2 - effMC1 * effMC2
                SF_den_up = effMC1_up + effMC2_up - effMC1_up * effMC2_up
                SF_den_down = effMC1_down + effMC2_down - effMC1_down * effMC2_down

                SF_num = effData1 + effData2 - effData1 * effData2
                SF_num_up = effData1_up + effData2_up - effData1_up * effData2_up
                SF_num_down = effData1_down + effData2_down - effData1_down * effData2_down

                SF = SF_num / SF_den
                SF_up = SF_num_up / SF_den_up
                SF_down = SF_num_down / SF_den_down

                sf_list.append(SF)
                sf_list_up.append(SF_up)
                sf_list_down.append(SF_down)

            HLT_weight = np.prod(sf_list) if sf_list else 1.0
            HLT_weightUp = np.prod(sf_list_up) if sf_list else 1.0
            HLT_weightDown = np.prod(sf_list_down) if sf_list else 1.0

            self.out.fillBranch("muTriggerWeight", HLT_weight)
            self.out.fillBranch("muTriggerWeightUp", HLT_weightUp)
            self.out.fillBranch("muTriggerWeightDown", HLT_weightDown)
        else:
            sf_list = []
            if len(event.selectedMuons) == 2:
                mu1, mu2 = event.selectedMuons[:2]

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
