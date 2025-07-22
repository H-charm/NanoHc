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

class MuonROOTProducer(Module):
    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]

        self.root_file = ROOT.TFile(f'../../data/MuonTRGSF/{self.year}/ScaleFactors_Muon_Z_HLT_{self.era}_eta_pt.root', "READ")
        # self.h_Mu_SF = self.root_file.Get("FINAL").Clone()
        self.h_Mu_SF = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt").Clone()   
        self.h_Mu_MCEff = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyMC").Clone()
        self.h_Mu_DataEff = self.root_file.Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt_efficiencyData").Clone()
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = self.dataset_type == "mc"
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('muHLTWeight', "F", limitedPrecision=10)
            self.out.branch('muTriggerWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muHLTWeightUp', "F", limitedPrecision=10)
                self.out.branch('muHLTWeightDown', "F", limitedPrecision=10)
                self.out.branch('muTriggerWeightUp', "F", limitedPrecision=10)
                self.out.branch('muTriggerWeightDown', "F", limitedPrecision=10)

            self.out.branch('muHLTWeightData', "F", limitedPrecision=10)
            self.out.branch('muTriggerWeightData', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muHLTWeightUpData', "F", limitedPrecision=10)
                self.out.branch('muHLTWeightDownData', "F", limitedPrecision=10)
                self.out.branch('muTriggerWeightUpData', "F", limitedPrecision=10)
                self.out.branch('muTriggerWeightDownData', "F", limitedPrecision=10)
    
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
        MCEff_list = []
        DataEff_list = []

        MCEff1_list = []
        DataEff1_list = []

        MCEff2_list = []
        DataEff2_list = []
            
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

            MCEff1_list.append(effMC1)
            MCEff2_list.append(effMC2)
            DataEff1_list.append(effData1)
            DataEff2_list.append(effData2)

            sf_list.append(SF)
            MCEff_list.append(SF_den)
            DataEff_list.append(SF_num)

        HLT_weight = np.prod(sf_list) if sf_list else 1.0

        HLT_MC_weight = np.prod(MCEff_list) if MCEff_list else 1.0
        HLT_Data_weight = np.prod(DataEff_list) if DataEff_list else 1.0

        self.out.fillBranch("muTriggerWeight", HLT_weight)

        self.out.fillBranch("muHLTWeight", HLT_MC_weight)
        self.out.fillBranch("muHLTWeightData", HLT_Data_weight)
        return True





# class ElectronSFProducer(Module):
#     def __init__(self, year, dataset_type, **kwargs):
#         self.year = year
#         self.dataset_type = dataset_type
#         self.era = era_dict[self.year]

#         base_path = f'../../data/ElectronSF/{self.year}/'
#         self.root_file_lowpt = ROOT.TFile(f'{base_path}egammaEffi_ptBelow20.root', "READ")
#         self.root_file_midpt = ROOT.TFile(f'{base_path}egammaEffi_ptBelow75.root', "READ")
#         self.root_file_highpt = ROOT.TFile(f'{base_path}egammaEffi_ptAbove75.root', "READ")
#         self.root_file_SF = ROOT.TFile(f'{base_path}SF2D_RMS.root', "READ")

#         self.h_el_lowpt = self.root_file_lowpt.Get("EGamma_SF2D").Clone()
#         self.h_el_midpt = self.root_file_midpt.Get("EGamma_SF2D").Clone()
#         self.h_el_highpt = self.root_file_highpt.Get("EGamma_SF2D").Clone()
#         self.h_el_SF = self.root_file_SF.Get("EGamma_SF2D").Clone()

#     def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#         self.isMC = self.dataset_type == "mc"
#         if self.isMC:
#             self.out = wrappedOutputTree
#             self.out.branch("elEffWeight", "F", limitedPrecision=12)

#     def analyze(self, event):
#         if not self.isMC:
#             return True

#         reco_sfs = []
#         id_sfs = []

#         for lep in event.selectedElectrons:
#             if abs(lep.pdgId) != 11:
#                 continue

#             sc_eta = lep.eta + getattr(lep, 'deltaEtaSC', 0.0)
#             sc_eta = max(min(sc_eta, 2.49), -2.49)
#             pt = min(lep.pt, 500.0)

#             if pt < 10:
#                 reco = 1.0
#                 id_ = self.getID_SF(pt, sc_eta)
#             else:
#                 reco = self.getRecoSF(pt, sc_eta)
#                 id_ = self.getID_SF(pt, sc_eta)

#             reco_sfs.append(reco)
#             id_sfs.append(id_)

#         reco_sf_prod = np.prod(reco_sfs) if reco_sfs else 1.0
#         id_sf_prod = np.prod(id_sfs) if id_sfs else 1.0
#         total_weight = reco_sf_prod * id_sf_prod

#         self.out.fillBranch("elEffWeight", total_weight)
#         return True

#     def getRecoSF(self, pt, eta):
#         if pt < 20:
#             h = self.h_el_lowpt
#         elif pt < 75:
#             h = self.h_el_midpt
#         else:
#             h = self.h_el_highpt

#         binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
#         binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
#         return h.GetBinContent(binX, binY)

#     def getID_SF(self, pt, eta):
#         h = self.h_el_SF
#         binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
#         binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
#         return h.GetBinContent(binX, binY)