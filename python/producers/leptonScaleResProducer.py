###
# Add branches for electron scale and resolution corrections.
# 
# Example: 2022EE MC smearing+uncertainties:
# eleSS = eleScaleRes("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS.json.gz", "Scale", "Smearing", True)
# 
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
from math import pi
import ROOT
from ROOT import MuonScaRe

era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}

class eleScaleRes(Module):
    def __init__(self, year, dataset_type, overwritePt=False, **kwargs):
        """Add branches for electron scale and resolution corrections.
        Parameters:
            json: full path of json file
            scaleKey: key for the scale correction for data and scale uncertainty in MC
            smearKey: key for the MC smearing correction (None for Data)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """
        self.year = year
        self.dataset_type = dataset_type
        self.overwritePt = overwritePt
        self.isMC = True if self.dataset_type == "mc" else False

        if self.year == "2022" or self.year == "2022EE":
            correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{era_dict[self.year]}/electronSS.json.gz'
            self.corr_scale = correctionlib.CorrectionSet.from_file(correction_file)["Scale"]
            self.corr_smear = correctionlib.CorrectionSet.from_file(correction_file)["Smearing"]
        elif self.year == "2023":
            correction_file = f"../../data/ElectronScale/electronSS_preBPix.json.gz"
            self.corr_scale = correctionlib.CorrectionSet.from_file(correction_file)["2023PromptC_ScaleJSON"]
            self.corr_smear = correctionlib.CorrectionSet.from_file(correction_file)["2023PromptC_SmearingJSON"]
        elif self.year == "2023BPix":
            correction_file = f"../../data/ElectronScale/electronSS_postBPix.json.gz"
            self.corr_scale = correctionlib.CorrectionSet.from_file(correction_file)["2023PromptD_ScaleJSON"]
            self.corr_smear = correctionlib.CorrectionSet.from_file(correction_file)["2023PromptD_SmearingJSON"]
        #evaluator = correctionlib.CorrectionSet.from_file(correction_file)
        
        

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Electron_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_uncorrected_pt", "F", lenVar="nElectron")
        else:
            self.out.branch("Electron_corrected_pt", "F", lenVar="nElectron")
        if self.isMC :
            self.out.branch("Electron_scaleUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_scaleDn_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearDn_pt", "F", lenVar="nElectron")

    def analyze(self, event):
        electrons = event.selectedElectrons

        pt_corr = []
        pt_smear_up = []
        pt_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []
        
        for ele in electrons:
            if self.isMC :
                # Set up a deterministic random seed.
                # The seed is unique by event and electron.
                # A fixed entropy value is also included to decorrelate different modules doing similar things.
                rng = np.random.default_rng(seed=np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((ele.phi/pi)%1)*1e12), 5402201385]))
                rho = self.corr_smear.evaluate("rho", ele.eta, ele.r9)
                smearing = rng.normal(loc=1., scale=rho)
                pt_corr.append(smearing * ele.pt)

                unc_rho = self.corr_smear.evaluate("err_rho", ele.eta, ele.r9)
                smearing_up = rng.normal(loc=1., scale=rho + unc_rho)
                smearing_dn = rng.normal(loc=1., scale=rho - unc_rho)
                pt_smear_up.append(smearing_up * ele.pt)
                pt_smear_dn.append(smearing_dn * ele.pt)

                scale_MC_unc = self.corr_scale.evaluate("total_uncertainty", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_scale_up.append((1+scale_MC_unc) * ele.pt)
                pt_scale_dn.append((1-scale_MC_unc) * ele.pt)
            else :
                scale = self.corr_scale.evaluate("total_correction", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_corr.append(scale * ele.pt)

        if self.overwritePt :
            pt_uncorr = list(ele.pt for ele in electrons)
            self.out.fillBranch("Electron_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Electron_pt", pt_corr)
        else :
            self.out.fillBranch("Electron_corrected_pt", pt_corr)

        if self.isMC :
            self.out.fillBranch("Electron_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Electron_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Electron_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Electron_scaleDn_pt", pt_scale_dn)

        return True

key_dict={
    "2022":     "2022_schemaV2.json.gz",
    "2022EE":   "2022EE_schemaV2.json.gz",
    "2023":     "2023_schemaV2.json.gz",
    "2023BPix": "2023BPix_schemaV2.json.gz"
}

class muonScaleRes(Module):
    def __init__(self, year, dataset_type, overwritePt=False, maxPt=200., **kwargs):
        """Add branches for muon scale and resolution corrections.
        Parameters:
            json: full path of json file
            isMC: True for MC (smear pt+add uncertainties), False for data (scale pt)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
            maxPt: # Do not correct for muons above this pT (200 GeV according to current recipe, 8/24)
        """

        self.year = year
        self.dataset_type = dataset_type
        self.overwritePt = overwritePt
        self.maxPt = maxPt
        self.isMC = True if self.dataset_type == "mc" else False
        
        correction_file = f"../../data/MuonScale/{key_dict[self.year]}"
        self.corrModule = MuonScaRe(correction_file)


    def getPtCorr(self, muon, var = "nom") :
        if muon.pt > self.maxPt : 
            return muon.pt
        isData = int(not self.isMC)
        scale_corr = self.corrModule.pt_scale(isData, muon.pt, muon.eta, muon.phi, muon.charge, var)
        pt_corr = scale_corr

        if self.isMC:
            smear_corr = self.corrModule.pt_resol(scale_corr, muon.eta, muon.nTrackerLayers, var)
            pt_corr = smear_corr

        return pt_corr


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Muon_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_uncorrected_pt", "F", lenVar="nMuon")
        else:
            self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon")
        if self.isMC:
            self.out.branch("Muon_syst_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_stat_pt", "F", lenVar="nMuon")


    def analyze(self, event):
        if event.selectedMuons == 0 :
            return True

        muons = event.selectedMuons 

        pt_corr = [0.]*len(muons)
        if self.isMC:
            pt_syst = [0.]*len(muons)
            pt_stat = [0.]*len(muons)
                
        for imu, muon in enumerate(muons):
            # Set up a deterministic random seed.
            # The seed is unique by event and muon.
            # A fixed entropy value is also included to decorrelate different modules doing similar things.
            seedSeq = np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((muon.phi/pi*100.)%1)*1e10), 351740215])
            self.corrModule.setSeed(int(seedSeq.generate_state(1,np.uint64)[0]))

            pt_corr[imu] = self.getPtCorr(muon, "nom")

            if self.isMC:
                # TODO: Check. Are we assuming up=dn?
                pt_syst[imu] = self.getPtCorr(muon, "syst")
                pt_stat[imu] = self.getPtCorr(muon, "stat")

        if self.overwritePt :
            pt_uncorr = list(mu.pt for mu in muons)
            self.out.fillBranch("Muon_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Muon_pt", pt_corr)
        else :
            self.out.fillBranch("Muon_corrected_pt", pt_corr)

        if self.isMC:
            self.out.fillBranch("Muon_syst_pt", pt_syst)
            self.out.fillBranch("Muon_stat_pt", pt_stat)

        return True