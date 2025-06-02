import numpy as np
import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

era_dict = {
    "2022": '2022_Summer22',
    "2022EE": '2022_Summer22EE',
    "2023": '2023_Summer23',
    "2023BPix": '2023_Summer23BPix'
}


class ElectronSFProducer(Module):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]

        base_path = f'../../data/ElectronSF/{self.year}/'
        self.root_file_lowpt = ROOT.TFile(f'{base_path}egammaEffi_ptBelow20.root', "READ")
        self.root_file_midpt = ROOT.TFile(f'{base_path}egammaEffi_ptBelow75.root', "READ")
        self.root_file_highpt = ROOT.TFile(f'{base_path}egammaEffi_ptAbove75.root', "READ")
        self.root_file_SF = ROOT.TFile(f'{base_path}SF2D_RMS.root', "READ")

        self.h_el_lowpt = self.root_file_lowpt.Get("EGamma_SF2D").Clone()
        self.h_el_midpt = self.root_file_midpt.Get("EGamma_SF2D").Clone()
        self.h_el_highpt = self.root_file_highpt.Get("EGamma_SF2D").Clone()
        self.h_el_SF = self.root_file_SF.Get("EGamma_SF2D").Clone()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = self.dataset_type == "mc"
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch("elEffWeight", "F", limitedPrecision=12)

    def analyze(self, event):
        if not self.isMC:
            return True

        reco_sfs = []
        id_sfs = []

        for lep in event.selectedElectrons:
            if abs(lep.pdgId) != 11:
                continue

            sc_eta = lep.eta + getattr(lep, 'deltaEtaSC', 0.0)
            sc_eta = max(min(sc_eta, 2.49), -2.49)
            pt = min(lep.pt, 500.0)

            if pt < 10:
                reco = 1.0
                id_ = self.getID_SF(pt, sc_eta)
            else:
                reco = self.getRecoSF(pt, sc_eta)
                id_ = self.getID_SF(pt, sc_eta)

            reco_sfs.append(reco)
            id_sfs.append(id_)

        reco_sf_prod = np.prod(reco_sfs) if reco_sfs else 1.0
        id_sf_prod = np.prod(id_sfs) if id_sfs else 1.0
        total_weight = reco_sf_prod * id_sf_prod

        self.out.fillBranch("elEffWeight", total_weight)
        return True

    def getRecoSF(self, pt, eta):
        if pt < 20:
            h = self.h_el_lowpt
        elif pt < 75:
            h = self.h_el_midpt
        else:
            h = self.h_el_highpt

        binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
        binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
        return h.GetBinContent(binX, binY)

    def getID_SF(self, pt, eta):
        h = self.h_el_SF
        binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
        binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
        return h.GetBinContent(binX, binY)


class MuonSFProducer(Module):
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]

        self.root_file = ROOT.TFile(f'../../data/MuonSF/{self.year}/Muon_corr.root', "READ")
        self.h_Mu_SF = self.root_file.Get("FINAL").Clone()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = self.dataset_type == "mc"
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch("muEffWeight", "F", limitedPrecision=12)

    def analyze(self, event):
        if not self.isMC:
            return True

        sf_list = []

        for lep in event.selectedMuons:
            if abs(lep.pdgId) != 13:
                continue

            pt = min(lep.pt, 199.0)
            eta = lep.eta
            sf = self.getSF(pt, eta)
            sf_list.append(sf)

        mu_eff_weight = np.prod(sf_list) if sf_list else 1.0

        self.out.fillBranch("muEffWeight", mu_eff_weight)
        return True

    def getSF(self, pt, eta):
        h = self.h_Mu_SF
        binX = min(max(1, h.GetXaxis().FindBin(eta)), h.GetNbinsX())
        binY = min(max(1, h.GetYaxis().FindBin(pt)), h.GetNbinsY())
        return h.GetBinContent(binX, binY)


