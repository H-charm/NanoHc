from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import ROOT
import math
import itertools
from functools import cmp_to_key
from ..helpers.utils import sumP4
ROOT.PyConfig.IgnoreCommandLineOptions = True

lumi_dict = {"2022": 7.9804, "2022EE": 26.6717, "2023": 17.794, "2023BPix": 9.451}

class Zcandidate:
    def __init__(self, lep1, lep2, fsrPhotons, fsrIndices):
        self.lep1 = lep1
        self.lep2 = lep2

        self.lep1_uncorrpt=lep1.pt
        self.lep2_uncorrpt=lep2.pt

        # Get fsr indices based on pdgId
        if abs(lep1.pdgId) == 13:
            self.fsrIdx1 = fsrIndices["muFsrPhotonIdx"][lep1.index]
            self.fsrIdx2 = fsrIndices["muFsrPhotonIdx"][lep2.index]
        elif abs(lep1.pdgId) == 11:
            self.fsrIdx1 = fsrIndices["eleFsrPhotonIdx"][lep1.index]
            self.fsrIdx2 = fsrIndices["eleFsrPhotonIdx"][lep2.index]
        else:
            self.fsrIdx1 = -1
            self.fsrIdx2 = -1

        # Get dressed four-momentum
        self.lep1_dressed = lep1.p4()
        self.lep2_dressed = lep2.p4()

        if self.fsrIdx1 >= 0:
            fsr1 = fsrPhotons[self.fsrIdx1]
            fsr1_p4 = ROOT.TLorentzVector()
            fsr1_p4.SetPtEtaPhiM(fsr1.pt, fsr1.eta, fsr1.phi, 0.0)
            self.lep1_dressed += fsr1_p4
        if self.fsrIdx2 >= 0:
            fsr2 = fsrPhotons[self.fsrIdx2]
            fsr2_p4 = ROOT.TLorentzVector()
            fsr2_p4.SetPtEtaPhiM(fsr2.pt, fsr2.eta, fsr2.phi, 0.0) 
            self.lep2_dressed += fsr2_p4

        self.p4 = self.lep1_dressed + self.lep2_dressed

        # Kinematic quantities from dressed p4
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()
        self.mass = self.p4.M()

    def sumpt(self):
        return self.lep1.pt + self.lep2.pt

    def finalState(self):
        return self.lep1.pdgId * self.lep2.pdgId

class ZZcandidate:

    def __init__(self,Z1,Z2):
        self.Z1 = Z1
        self.Z2 = Z2
        self.pt = sumP4(self.Z1, self.Z2).Pt()
        self.eta = sumP4(self.Z1, self.Z2).Eta()
        self.phi = sumP4(self.Z1, self.Z2).Phi()
        self.mass = sumP4(self.Z1, self.Z2).M()

class ZLLcandidate:

    def __init__(self,Z1,Z2):
        self.Z1 = Z1
        self.Z2 = Z2
        self.pt = sumP4(self.Z1, self.Z2).Pt()
        self.eta = sumP4(self.Z1, self.Z2).Eta()
        self.phi = sumP4(self.Z1, self.Z2).Phi()
        self.mass = sumP4(self.Z1, self.Z2).M()
        self.category = None 

        self.lep3_pt = Z2.lep1.pt
        self.lep3_eta = Z2.lep1.eta
        self.lep3_pdgId = Z2.lep1.pdgId

        self.lep4_pt = Z2.lep2.pt
        self.lep4_eta = Z2.lep2.eta
        self.lep4_pdgId = Z2.lep2.pdgId 

        self.lep3_lostHits = Z2.lep1.lostHits if abs(Z2.lep1.pdgId)==11 else -1
        self.lep4_lostHits = Z2.lep2.lostHits if abs(Z2.lep2.pdgId)==11 else -1

        self.lep3_isGood= Z2.lep1.isFullID
        self.lep4_isGood= Z2.lep2.isFullID


class ZLcandidate:

    def __init__(self,Z1,lep):
        self.Z1 = Z1
        self.lep = lep
        self.pt = sumP4(self.Z1, self.lep).Pt()
        self.eta = sumP4(self.Z1, self.lep).Eta()
        self.phi = sumP4(self.Z1, self.lep).Phi()
        p1=ROOT.TLorentzVector()
        p2=ROOT.TLorentzVector()
        p3=ROOT.TLorentzVector()
        p1.SetPtEtaPhiM(Z1.lep1.pt, Z1.lep1.eta, Z1.lep1.phi, 0.0)
        p2.SetPtEtaPhiM(Z1.lep2.pt, Z1.lep2.eta, Z1.lep2.phi, 0.0)
        p3.SetPtEtaPhiM(lep.pt, lep.eta, lep.phi, 0.0)
        self.trimass = ((p1+p2)+p3).M()
        self.mass = self.Z1.mass
        self.pt2 = lep.pt
        self.pdgId = lep.pdgId
        self.lostHits = lep.lostHits if abs(lep.pdgId)==11 else -1
        self.phi2= lep.phi
        self.eta2 = lep.eta
        self.sip3d=lep.sip3d
        self.dxy=lep.dxy
        self.dz=lep.dz
        self.pfcand=lep.isPFcand


class CRProducer(Module):
    
    def __init__(self, year, dataset_type, sample):
        self.year = year
        self.sample = sample
        self.dataset_type = dataset_type

        # Define the variables you want to plot
        self.lep_vars = ["pt","eta","phi","pdgId"]
        self.mu_vars = ["pt","eta","phi", "pdgId"]
        self.jet_vars = ["pt","eta","phi","mass","bdisc","cvbdisc","cvldisc","gvudsdisc"]
        self.jet_vars_mc = ["hadronFlavour"]
        self.Z_vars = ["pt","eta","phi","mass","onshell_mass","offshell_mass"]
        self.ZZ_vars = ["pt","eta","phi","mass"] 
        self.H_vars = ["pt","eta","phi","mass"]  
        self.H4e_vars=["pt","eta","phi","mass"]
        self.H4mu_vars=["pt","eta","phi","mass"]
        self.H2e2mu_vars=["pt","eta","phi","mass"]  
        self.ZZ4e_vars=["pt","eta","phi","mass"]
        self.ZZ4mu_vars=["pt","eta","phi","mass"]
        self.ZZ2e2mu_vars=["pt","eta","phi","mass"]

        self.ZLL_vars=["pt","eta","phi","mass","lep3_pt","lep3_eta","lep3_pdgId","lep4_pt","lep4_eta","lep4_pdgId","lep3_lostHits","lep4_lostHits","lep3_isGood","lep4_isGood"] 
        self.ZL_vars=["pt","eta","phi","mass","trimass","pt2", "eta2","sip3d","dz","dxy","iso","pfcand","pdgId","lostHits","phi2"]

        # Define the prefixes
        self.mu_prefix = "mu_"
        self.el_prefix = "el_"
        self.lep_prefix = "lep_"
        self.full_mu_prefix = "full_mu_"
        self.full_el_prefix = "full_el_"
        self.full_lep_prefix = "full_lep_"
        self.jet_prefix = "jet_"
        self.Z_prefix = "Z_"
        self.ZZ_prefix = "ZZ_"
        self.ZZ4e_prefix = "ZZ4e_"
        self.ZZ4mu_prefix = "ZZ4mu_"
        self.ZZ2e2mu_prefix = "ZZ2e2mu_"
        self.H_prefix = "H_"
        self.H4e_prefix = "H4e_"
        self.H4mu_prefix = "H4mu_"
        self.H2e2mu_prefix = "H2e2mu_"

        self.ZLall_prefix = "ZLall_"
        self.ZLalle_prefix = "ZLalle_"
        self.ZLallmu_prefix = "ZLallmu_"
        self.ZLpass_prefix = "ZLpass_"
        self.ZLpasse_prefix = "ZLpasse_"
        self.ZLpassmu_prefix = "ZLpassmu_"

        self.ZLL_prefix = "ZLL_"
        self.ZLL4e_prefix = "ZLL4e_"
        self.ZLL4mu_prefix = "ZLL4mu_"
        self.ZLL2e2mu_prefix = "ZLL2e2mu_"
        self.ZLL2mu2e_prefix = "ZLL2mu2e_"

        self.ZLL2P2F_prefix = "ZLL2P2F_"
        self.ZLL2P2F4e_prefix = "ZLL2P2F4e_"
        self.ZLL2P2F4mu_prefix = "ZLL2P2F4mu_"
        self.ZLL2P2F2e2mu_prefix = "ZLL2P2F2e2mu_"
        self.ZLL2P2F2mu2e_prefix = "ZLL2P2F2mu2e_"

        self.ZLL3P1F_prefix = "ZLL3P1F_"
        self.ZLL3P1F4e_prefix = "ZLL3P1F4e_"
        self.ZLL3P1F4mu_prefix = "ZLL3P1F4mu_"
        self.ZLL3P1F2e2mu_prefix = "ZLL3P1F2e2mu_"
        self.ZLL3P1F2mu2e_prefix = "ZLL3P1F2mu2e_"

        self.ZLLSSCR_prefix = "ZLLSSCR_"
        self.ZLLSSCR4e_prefix = "ZLLSSCR4e_"
        self.ZLLSSCR4mu_prefix = "ZLLSSCR4mu_"
        self.ZLLSSCR2e2mu_prefix = "ZLLSSCR2e2mu_"
        self.ZLLSSCR2mu2e_prefix = "ZLLSSCR2mu2e_"

        self.ZZSR_prefix = "ZZSR_"
        self.ZZSR4e_prefix = "ZZSR4e_"
        self.ZZSR4mu_prefix = "ZZSR4mu_"
        self.ZZSR2e2mu_prefix = "ZZSR2e2mu_"

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        
        self.out = wrappedOutputTree

        for lep_var in self.lep_vars:
            self.out.branch(self.mu_prefix + lep_var, "F", 20, lenVar="nMu")
            self.out.branch(self.el_prefix + lep_var, "F", 20, lenVar="nEl")
            self.out.branch(self.lep_prefix + lep_var, "F", 20, lenVar="nLep")
            self.out.branch(self.full_mu_prefix + lep_var, "F", 20, lenVar="nfullMu")
            self.out.branch(self.full_el_prefix + lep_var, "F", 20, lenVar="nfullEl")
            self.out.branch(self.full_lep_prefix + lep_var, "F", 20, lenVar="nfullLep")
        
        # Define jet branches
        for jet_var in self.jet_vars:
            self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar="nJet")
        if self.isMC: 
            for jet_var in self.jet_vars_mc:
                self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar="nJet")
        
        # Define trigger branches
        self.out.branch("HLT_passZZ4lEle", "O")   # pass Ele triggers
        self.out.branch("HLT_passZZ4lMu", "O")    # pass Muon triggers
        self.out.branch("HLT_passZZ4lMuEle", "O") # pass MuEle triggers
        self.out.branch("HLT_passZZ4l", "O")      # pass trigger requirements for the given sample (including sample precedence vetos) 

        # Define luminosity branch
        self.out.branch("lumiwgt", "F")

        # Define branches for the Zcandidates
        for Z_var in self.Z_vars:
            self.out.branch(self.Z_prefix + Z_var, "F", 20, lenVar="nZ")

        # Define branches for the ZZcandidates
        for ZZ_var in self.ZZ_vars:
            self.out.branch(self.ZZ_prefix + ZZ_var, "F", 20, lenVar="nZZ")

        for ZZ4e_var in self.ZZ4e_vars:
            self.out.branch(self.ZZ4e_prefix + ZZ4e_var, "F", 20, lenVar="nZZ4e")

        for ZZ4mu_var in self.ZZ4mu_vars:
            self.out.branch(self.ZZ4mu_prefix + ZZ4mu_var, "F", 20, lenVar="nZZ4mu")

        for ZZ2e2mu_var in self.ZZ2e2mu_vars:
            self.out.branch(self.ZZ2e2mu_prefix + ZZ2e2mu_var, "F", 20, lenVar="nZZ2e2mu")

        # Define branches for the ZLLcandidates
        for ZLL_var in self.ZLL_vars:
            self.out.branch(self.ZLL_prefix + ZLL_var, "F", 20, lenVar="nZLL" )
            self.out.branch(self.ZLL4e_prefix + ZLL_var, "F", 20, lenVar="nZLL4e" )
            self.out.branch(self.ZLL4mu_prefix + ZLL_var, "F", 20, lenVar="nZLL4mu" )
            self.out.branch(self.ZLL2e2mu_prefix + ZLL_var, "F", 20, lenVar="nZLL2e2mu" )
            
            self.out.branch(self.ZLL2P2F_prefix + ZLL_var, "F", 20, lenVar="nZLL2P2F" ) 
            self.out.branch(self.ZLL2P2F4e_prefix + ZLL_var, "F", 20, lenVar="nZLL2P2F4e" )
            self.out.branch(self.ZLL2P2F4mu_prefix + ZLL_var, "F", 20, lenVar="nZLL2P2F4mu" )
            self.out.branch(self.ZLL2P2F2e2mu_prefix + ZLL_var, "F", 20, lenVar="nZLL2P2F2e2mu" )

            self.out.branch(self.ZLL3P1F_prefix + ZLL_var, "F", 20, lenVar="nZLL3P1F" )
            self.out.branch(self.ZLL3P1F4e_prefix + ZLL_var, "F", 20, lenVar="nZLL3P1F4e" )
            self.out.branch(self.ZLL3P1F4mu_prefix + ZLL_var, "F", 20, lenVar="nZLL3P1F4mu" )
            self.out.branch(self.ZLL3P1F2e2mu_prefix + ZLL_var, "F", 20, lenVar="nZLL3P1F2e2mu" )

            self.out.branch(self.ZLLSSCR_prefix + ZLL_var, "F", 20, lenVar="nZLLSSCR" )
            self.out.branch(self.ZLLSSCR4e_prefix + ZLL_var, "F", 20, lenVar="nZLLSSCR4e" )
            self.out.branch(self.ZLLSSCR4mu_prefix + ZLL_var, "F", 20, lenVar="nZLLSSCR4mu" )
            self.out.branch(self.ZLLSSCR2e2mu_prefix + ZLL_var, "F", 20, lenVar="nZLLSSCR2e2mu" )

        for ZL_var in self.ZL_vars:
            self.out.branch(self.ZLall_prefix + ZL_var, "F", 20, lenVar="nZLall" )
            self.out.branch(self.ZLalle_prefix + ZL_var, "F", 20, lenVar="nZLalle" )
            self.out.branch(self.ZLallmu_prefix + ZL_var, "F", 20, lenVar="nZLallmu" )
            self.out.branch(self.ZLpass_prefix + ZL_var, "F", 20, lenVar="nZLpass" )
            self.out.branch(self.ZLpasse_prefix + ZL_var, "F", 20, lenVar="nZLpasse" )
            self.out.branch(self.ZLpassmu_prefix + ZL_var, "F", 20, lenVar="nZLpassmu" )
        # self.out.branch(self.ZL_prefix + "index", "I", 20, lenVar="nZL")

        for ZZ_var in self.ZZ_vars:
            self.out.branch(self.ZZSR_prefix + ZZ_var, "F", 20, lenVar="nZZSR")
            self.out.branch(self.ZZSR4e_prefix + ZZ_var, "F", 20, lenVar="nZZSR4e")
            self.out.branch(self.ZZSR4mu_prefix + ZZ_var, "F", 20, lenVar="nZZSR4mu")
            self.out.branch(self.ZZSR2e2mu_prefix + ZZ_var, "F", 20, lenVar="nZZSR2e2mu")

        # Define branches for the Hcandidates
        for H_var in self.H_vars:
            self.out.branch(self.H_prefix + H_var, "F", 20, lenVar="nH")

        for H4e_var in self.H4e_vars:
            self.out.branch(self.H4e_prefix + H4e_var, "F", 20, lenVar="nH4e")

        for H4mu_var in self.H4mu_vars:
            self.out.branch(self.H4mu_prefix + H4mu_var, "F", 20, lenVar="nH4mu")

        for H2e2mu_var in self.H2e2mu_vars:
            self.out.branch(self.H2e2mu_prefix + H2e2mu_var, "F", 20, lenVar="nH2e2mu")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if event.PV_npvsGood < 1: return False

        if event.MET_pt > 25: return False

        # Apply trigger selections 
        if self._select_triggers(event) is False:
            return False
        self._associate_fsr_photons_and_leptons(event)
        event.fullIDLeptons = event.fullIDMuons + event.fullIDElectrons
        # if len(event.fullIDLeptons) < 4: return False 
        event.selectedLeptons = event.relaxedMuons + event.relaxedElectrons

        if len(event.selectedLeptons) < 3: return False # For Z+L CR we need to change this to 3 

        self._select_jets(event)
        
        self._select_Z_candidates(event)

        self._build_SR_and_CR_combinations(event)
        
        # self._select_ZZ_candidates(event)
        # self._select_H_candidates(event)
        
        self._fill_event_info(event)

        return True
        
    def _select_triggers(self, event):

        passTrigger = False 
        out_data = {}
        if self.year == "2022" or self.year == "2022EE" : # Checked that these are unprescaled in run 359751
            passSingleEle = event.HLT_Ele30_WPTight_Gsf #Note: we used Ele32 in 2018! 
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        elif self.year == "2023" or self.year == "2023BPix" : # Checked that these are unprescaled, reference twikis for 2023 Eg & Muon Triggers https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary & https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2023
            passSingleEle = event.HLT_Ele30_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        else:
            print(f"Year {self.year} not found")

        if self.isMC or self.sample == "any" :
            passTrigger = passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu
        else: # Data: ensure each event is taken only from a single sample
            if self.sample == "" : sys.exit("ERROR: sample must be set in data") # we may want to merge triggers for test runs 
            if (self.sample in ["DoubleEle", "DoubleEG", "EGamma"] and (passDiEle or passTriEle)) or \
            (self.sample in ["Muon", "DoubleMu", "DoubleMuon"] and (passDiMu or passTriMu) and not passDiEle and not passTriEle) or \
            (self.sample in ["MuEG", "MuonEG"] and passMuEle and not passDiEle and not passTriEle and not passDiMu and not passTriMu) or \
            (self.sample in ["SingleElectron", "EGamma"] and passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
            (self.sample in ["SingleMuon", "Muon"] and passSingleMu and not passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle):
                passTrigger = True

        self.out.fillBranch("HLT_passZZ4lEle", passSingleEle or passDiEle or passTriEle)
        self.out.fillBranch("HLT_passZZ4lMu", passSingleMu or passDiMu or passTriMu)
        self.out.fillBranch("HLT_passZZ4lMuEle", passMuEle)
        self.out.fillBranch("HLT_passZZ4l", passTrigger)

        return passTrigger

    def _select_Z_candidates(self, event):

        event.Zcandidates = []
        event.ZcandidatesSR =[]
        event.bestZIdx = -1
        bestMassDiff = 9999

        ZmassNominal = 91.1876  # Z mass in GeV
        lepton_pairs = list(itertools.combinations(event.selectedLeptons, 2))

        for lep1, lep2 in lepton_pairs:
            if not abs(lep1.pdgId)==abs(lep2.pdgId):
                continue

            # Ensure consistent ordering (lep1 = negative charge)
            lep1, lep2 = (lep1, lep2) if lep1.pdgId < 0 else (lep2, lep1)

            fsrPhotons = Collection(event, "FsrPhoton")
            fsrIndices = {
                "muFsrPhotonIdx": getattr(event, "muFsrPhotonIdx", []),
                "eleFsrPhotonIdx": getattr(event, "eleFsrPhotonIdx", [])
            }
            # if not ((lep1.pt > 20 or lep2.pt > 20) and (lep1.pt > 10 or lep2.pt > 10)):
            #     continue

            Zcand = Zcandidate(lep1, lep2, fsrPhotons, fsrIndices)

            num_passed_pt20 = 0
            num_passed_pt10 = 0

            for lep in [lep1,lep2]:
                if lep.pt > 20:
                        num_passed_pt20 += 1
                elif lep.pt > 10:
                        num_passed_pt10 += 1

            if num_passed_pt20 == 0 or (num_passed_pt10 + num_passed_pt20) < 2:
                continue

            if Zcand.mass < 12 or Zcand.mass > 120:
                continue


            # Initialize flags
            Zcand.isOSSF = False
            Zcand.isSR = False
            Zcand.is1FCR = False
            Zcand.is2FCR = False
            Zcand.isSSCR = False

            # Lepton ID checks
            relaxed1 = lep1.isRelaxed
            relaxed2 = lep2.isRelaxed
            full1 = lep1.isFullID
            full2 = lep2.isFullID

            if (lep1.pdgId + lep2.pdgId) == 0:
                Zcand.isOSSF = True
                if relaxed1 and relaxed2:
                    nFull = int(full1) + int(full2)
                    if nFull == 2:
                        Zcand.isSR = True
                    elif nFull == 1:
                        Zcand.is1FCR = True
                        if full1 and not full2:
                            Zcand.leptons = (lep2, lep1)
                        elif not full1 and full2:
                            Zcand.leptons = (lep1, lep2)
                    elif nFull == 0:
                        Zcand.is2FCR = True

            # Same-sign control region
            elif lep1.pdgId == lep2.pdgId:
                if relaxed1 and relaxed2:
                    Zcand.isSSCR = True

            elif abs(lep1.pdgId) != abs(lep2.pdgId):
                continue

            # Save Z candidate
            event.Zcandidates.append(Zcand) 

            if Zcand.isSR:
                event.ZcandidatesSR.append(Zcand) 

            # Track best SR Z candidate
            if Zcand.isSR:
                massDiff = abs(Zcand.mass - ZmassNominal)
                if massDiff < bestMassDiff:
                    event.bestZIdx = len(event.Zcandidates) - 1
                    bestMassDiff = massDiff


    def _build_SR_and_CR_combinations(self, event):

        # Initialize storage
        event.ZZSRcandidates, event.ZLLcandidates = [], []
        event.ZLcandidates = None
        event.ZLcandidates_all, event.ZLcandidates_alle, event.ZLcandidates_allmu = [], [], []
        event.ZLcandidates_pass, event.ZLcandidates_passe, event.ZLcandidates_passmu = [], [], []
        event.bestZZIdx = -1

        Zs = event.Zcandidates
        ZZs, ZLLs, ZLLsTemp = [], [], []
        ZLL2P2F=[]
        ZLL3P1F=[]
        ZLLSSCR=[]
        ZLs_all = []
        ZLs_alle = []
        ZLs_allmu = []
        ZLs_pass = []
        ZLs_passe = []
        ZLs_passmu = []
        event.ZLL2P2Fcandidates = []
        event.ZLL3P1Fcandidates = []
        event.ZLLSSCRcandidates = []

        bestZZIdx, best2P2FCRIdx, best3P1FCRIdx, bestSSCRIdx = -1, -1, -1, -1
        ZmassNominal = 91.1876

        def bestCandCmp(a, b):
            if abs(a.Z1.mass - b.Z1.mass) < 1e-4:
                return -1 if a.Z2.lep1.pt + a.Z2.lep2.pt > b.Z2.lep1.pt + b.Z2.lep2.pt else 1
            return -1 if abs(a.Z1.mass - ZmassNominal) < abs(b.Z1.mass - ZmassNominal) else 1

        def leptonOverlapOrFail(leps):
            if len(set(leps)) != 4:
                return True
            lep_pts = sorted([lep.pt for lep in leps])
            if not (lep_pts[3] > 20 and lep_pts[2] > 10):
                return True
            lepton_pairs = itertools.combinations(leps, 2)
            dr_ll = [deltaR(l1.eta, l1.phi, l2.eta, l2.phi) for l1, l2 in lepton_pairs]
            if any(dr <= 0.02 for dr in dr_ll):
                return True
            m_ll = [(l1.p4() + l2.p4()).M() for l1, l2 in lepton_pairs if l1.pdgId + l2.pdgId == 0]
            return any(m <= 4 for m in m_ll)

        def failsSmartPairing(Z1, Z2):
            if abs(Z1.lep1.pdgId) != abs(Z2.lep1.pdgId): return False
            fsrPhotons = Collection(event, "FsrPhoton")
            fsrIndices = {
                "muFsrPhotonIdx": getattr(event, "muFsrPhotonIdx", []),
                "eleFsrPhotonIdx": getattr(event, "eleFsrPhotonIdx", [])
            }
            Zt1 = Zcandidate(Z1.lep1, Z2.lep2, fsrPhotons, fsrIndices)
            Zt2 = Zcandidate(Z2.lep1, Z1.lep2, fsrPhotons, fsrIndices)
            Za, Zb = (Zt1, Zt2) if abs(Zt1.mass - ZmassNominal) < abs(Zt2.mass - ZmassNominal) else (Zt2, Zt1)
            return abs(Za.mass - ZmassNominal) < abs(Z1.mass - ZmassNominal) and Zb.mass < 12

        # ---------- Signal Region (ZZ) ----------
        for i, Z1 in enumerate(Zs):
            for j in range(i + 1, len(Zs)):
                Z2 = Zs[j]
                if not (Z1.isOSSF and Z2.isOSSF): continue
                if any(l in [Z1.lep1, Z1.lep2] for l in [Z2.lep1, Z2.lep2]): continue
                if leptonOverlapOrFail([Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2]): continue
                if failsSmartPairing(Z1, Z2): continue

                ZZ = ZZcandidate(Z1, Z2)
                if ZZ.mass < 70: continue

                if Z1.isSR and Z2.isSR:
                    ZZ.isSR = True
                    ZZs.append(ZZ)
                    if bestZZIdx < 0 or bestCandCmp(ZZ, ZZs[bestZZIdx]) < 0:
                        bestZZIdx = len(ZZs) - 1

        event.ZZSRcandidates = ZZs
        if bestZZIdx >= 0:
            event.bestZZIdx = bestZZIdx

        # ---------- Control Regions (ZLL) ----------
        if bestZZIdx < 0 and len(Zs) >= 2:
            for Z1 in Zs:
                if not Z1.isSR: continue
                for Z2 in Zs:
                    if Z1 == Z2 or not (Z2.is1FCR or Z2.is2FCR or Z2.isSSCR): continue
                    if any(l in [Z1.lep1, Z1.lep2] for l in [Z2.lep1, Z2.lep2]): continue
                    if Z1.mass < 40: continue
                    if leptonOverlapOrFail([Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2]): continue
                    if failsSmartPairing(Z1, Z2): continue

                    ZLL = ZLLcandidate(Z1, Z2)
                    # if ZLL.mass < 100: continue
                    ZLL.Z2 = Z2
                    ZLLsTemp.append(ZLL)

                    if Z2.is2FCR and (best2P2FCRIdx < 0 or bestCandCmp(ZLL, ZLLsTemp[best2P2FCRIdx]) < 0):
                        best2P2FCRIdx = len(ZLLsTemp) - 1
                        ZLL.category = "2P2F"
                    if Z2.is1FCR and (best3P1FCRIdx < 0 or bestCandCmp(ZLL, ZLLsTemp[best3P1FCRIdx]) < 0):
                        best3P1FCRIdx = len(ZLLsTemp) - 1
                        ZLL.category = "3P1F"
                    if Z2.isSSCR and (bestSSCRIdx < 0 or bestCandCmp(ZLL, ZLLsTemp[bestSSCRIdx]) < 0):
                        bestSSCRIdx = len(ZLLsTemp) - 1
                        ZLL.category = "SSCR"

            # Keep only best ZLL candidates per category
            for idx, ZLL in enumerate(ZLLsTemp):
                if idx in [best2P2FCRIdx, best3P1FCRIdx, bestSSCRIdx]:
                    ZLLs.append(ZLL)
                if idx == best2P2FCRIdx:
                    ZLL2P2F.append(ZLL)
                if idx == best3P1FCRIdx:
                    ZLL3P1F.append(ZLL)
                if idx == bestSSCRIdx:
                    ZLLSSCR.append(ZLL)

        event.ZLLcandidates = ZLLs
        event.ZLL2P2Fcandidates = ZLL2P2F
        event.ZLL3P1Fcandidates = ZLL3P1F
        event.ZLLSSCRcandidates = ZLLSSCR

         # ---------- Control Region: Z + L ----------
        if hasattr(event, "bestZIdx") and event.bestZIdx >= 0:
            bestZ = event.Zcandidates[event.bestZIdx]
            if (bestZ.mass > 40) and (bestZ.mass < 120): 
                # extra_leptons = [
                #     lep for lep in event.selectedLeptons
                #     if lep not in [bestZ.lep1, bestZ.lep2] and lep.isRelaxed
                # ]
                selected_with_idx = list(enumerate(event.selectedLeptons))

                # Find index of bestZ.lep1 and bestZ.lep2
                idx_lep1 = next((i for i, lep in selected_with_idx if lep is bestZ.lep1), -1)
                idx_lep2 = next((i for i, lep in selected_with_idx if lep is bestZ.lep2), -1)

                # Filter out Z leptons and keep relaxed leptons
                extra_leptons_with_idx = [
                    (i, lep) for i, lep in selected_with_idx
                    if lep not in [bestZ.lep1, bestZ.lep2] and lep.isRelaxed
                ]
                if len(extra_leptons_with_idx) == 1:
                    idx_aL, aL = extra_leptons_with_idx[0]
                    # print("DEBUG: Z+L candidate leptons:")
                    # print(f" - aL (third lepton): index={idx_aL}, pt={aL.pt}")

                    if deltaR(aL.eta, aL.phi, bestZ.lep1.eta, bestZ.lep1.phi) > 0.02 and \
                        deltaR(aL.eta, aL.phi, bestZ.lep2.eta, bestZ.lep2.phi) > 0.02:

                        # QCD suppression (m_ll > 4 GeV if OS)
                        if (aL.charge != bestZ.lep1.charge and sumP4(aL, bestZ.lep1).M() <= 4) or \
                            (aL.charge != bestZ.lep2.charge and sumP4(aL, bestZ.lep2).M() <= 4):
                            pass
                        else:
                            event.ZLcandidates = aL
                            ZL = ZLcandidate(bestZ, aL)                 
                            ZLs_all.append(ZL)
                            if abs(aL.pdgId) == 11:
                                ZLs_alle.append(ZL)
                            elif abs(aL.pdgId) == 13:
                                ZLs_allmu.append(ZL)

                            if aL.isFullID:
                                ZLs_pass.append(ZL)
                                if abs(aL.pdgId) == 11:
                                    ZLs_passe.append(ZL)
                                elif abs(aL.pdgId) == 13:
                                    ZLs_passmu.append(ZL)
        
        event.ZLcandidates_all = ZLs_all
        event.ZLcandidates_alle = ZLs_alle
        event.ZLcandidates_allmu = ZLs_allmu
        event.ZLcandidates_pass = ZLs_pass
        event.ZLcandidates_passe = ZLs_passe
        event.ZLcandidates_passmu = ZLs_passmu



    def _flag_onshell_and_offshell_Z(self, Zcand_pair):
                    
        Z1 = Zcand_pair[0]
        Z2 = Zcand_pair[1]
        
        mZ = 91.1876
        
        ## minimal |mZcandidate - mZ| is considered on-shell
        d_mZ1_mZ = abs( Z1.mass - mZ )
        d_mZ2_mZ = abs( Z2.mass - mZ )
        onshell_idx = 0 if d_mZ1_mZ < d_mZ2_mZ else 1
        offshell_idx = 1 if d_mZ1_mZ < d_mZ2_mZ else 0
        
        Zcand_pair[onshell_idx].is_onshell = True
        Zcand_pair[offshell_idx].is_onshell = False
                
    def _select_ZZ_candidates(self, event):

        event.ZZcandidates = []

        if  len(event.ZcandidatesSR) < 2: return False 
                
        Zcand_pairs = list(itertools.combinations(event.ZcandidatesSR, 2))
        for Zcand_pair in Zcand_pairs:
                        
            self._flag_onshell_and_offshell_Z(Zcand_pair) 
            Z1, Z2 = (Zcand_pair[0], Zcand_pair[1]) if Zcand_pair[0].is_onshell else (Zcand_pair[1], Zcand_pair[0]) # by definition Z1 onshell, Z2 offshell

            lep_ids = {
                abs(Z1.lep1.pdgId),
                abs(Z1.lep2.pdgId),
                abs(Z2.lep1.pdgId),
                abs(Z2.lep2.pdgId)
            }

   
            if Z1.mass < 40: continue  
            # if any(l in [Z1.lep1, Z1.lep2] for l in [Z2.lep1, Z2.lep2]):
            #     continue
            
            leptons = [Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2] 
            
            ## two DISTINCT leptons must pass pt > 10 and pt > 20
            # at_least_one_passed_pt20 = False
            # at_least_one_passed_pt10 = False
            # for lep in leptons:
            #     if lep.pt > 20: at_least_one_passed_pt20 = True
            #     elif lep.pt > 10: at_least_one_passed_pt10 = True
            # if not (at_least_one_passed_pt10 and at_least_one_passed_pt20): continue

            # Two DISTINCT leptons must pass pt > 10 and pt > 20
            # num_passed_pt20 = 0
            # num_passed_pt10 = 0

            # for lep in leptons:
            #     if lep.pt > 20:
            #         num_passed_pt20 += 1
            #     elif lep.pt > 10:
            #         num_passed_pt10 += 1

            # # Ensure we have at least one lepton passing each threshold, and they are distinct
            # if num_passed_pt20 == 0 or (num_passed_pt10 + num_passed_pt20) < 2: 
            #     print("pt cut!")
            #     continue 
            lep_pts = sorted([lep.pt for lep in leptons])
            if not (lep_pts[3] > 20 and lep_pts[2] > 10):
                continue
            
            lepton_pairs = list(itertools.combinations(leptons, 2)) # 6 combinations
            
            dr_ll_values = [] # between each of the four leptons (Ghost removal: all pairs)
            m_ll_values = [] # between each of the four leptons (QCD suppression: opposite-sign pairs and same flavour)                     
            for lepton_pair in lepton_pairs:
                lep1 = lepton_pair[0]
                lep2 = lepton_pair[1]
                            
                dr = deltaR(lep1.eta,lep1.phi,lep2.eta,lep2.phi)
                dr_ll_values.append(dr)
                
                if (lep1.pdgId + lep2.pdgId) == 0:
                    m_ll = (lep1.p4()+lep2.p4()).M()
                    m_ll_values.append(m_ll)
                    
            if any(dr <= 0.02 for dr in dr_ll_values): 
                continue
            if any(m_ll <= 4 for m_ll in m_ll_values): 
                continue
            
            ## smart cut: check alternative pairing (4e or 4mu)
            if abs(Z1.lep1.pdgId) == abs(Z1.lep2.pdgId) == abs(Z2.lep1.pdgId) == abs(Z2.lep2.pdgId):

                fsrPhotons = Collection(event, "FsrPhoton")
                fsrIndices = {
                    "muFsrPhotonIdx": getattr(event, "muFsrPhotonIdx", []),
                    "eleFsrPhotonIdx": getattr(event, "eleFsrPhotonIdx", [])
                }

                Ztemp1 = Zcandidate(Z1.lep1, Z2.lep2, fsrPhotons, fsrIndices)
                Ztemp2 = Zcandidate(Z2.lep1, Z1.lep2, fsrPhotons, fsrIndices)
                mZ = 91.1876
                Za, Zb = (Ztemp1, Ztemp2) if ( abs(Ztemp1.mass - mZ) < abs(Ztemp2.mass - mZ) ) else (Ztemp2, Ztemp1)
                                
                ## reject if |mZa - mZ| < |mZ1 - mZ| and mZb < 12
                if ( abs(Za.mass - mZ) < abs(Z1.mass - mZ) ) and Zb.mass < 12: 
                    continue
                    
            ZZcand = ZZcandidate(Z1,Z2)
            m_4l = ZZcand.mass
            if m_4l < 70: continue
            event.ZZcandidates.append(ZZcand)      
            
    def _select_H_candidates(self, event):
        
        def best_candidate_comparator(a, b):
            if abs(a.Z1.mass - b.Z1.mass) < 1e-4 : # If Z1 masses are similar, compare Z2 sum of transverse momenta
                if a.Z2.lep1.pt + a.Z2.lep2.pt > b.Z2.lep1.pt + b.Z2.lep2.pt:
                    return -1  # a is better
                else:
                    return 1   # b is better
            else:
                if abs(a.Z1.mass - 91.1876) < abs(b.Z1.mass - 91.1876):
                    return -1  # a is better
                else:
                    return 1   # b is better days ago

        event.Hcandidates = []

        if len(event.ZZcandidates) == 0:
            return
        ##---------------------
        ## Comparator
        best_candidate = min(event.ZZcandidates, key=cmp_to_key(best_candidate_comparator))
        event.Hcandidates.append(best_candidate) 

    def isoFsrCorr(self, lep, selectedFSR):
        combRelIsoPFFSRCorr = lep.pfRelIso03_all
        for fsr in selectedFSR:
            dR = deltaR(lep.eta, lep.phi, fsr.eta, fsr.phi)
            if 0.01 < dR < 0.3:
                combRelIsoPFFSRCorr = max(0.0, combRelIsoPFFSRCorr - fsr.pt / lep.pt)
        return combRelIsoPFFSRCorr

        
    ## taken from here https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/nanoZZ4lAnalysis.py    
    def _associate_fsr_photons_and_leptons(self, event):
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        fsrPhotons = Collection(event, "FsrPhoton")

        fsr_dRET2_max = 0.012
        relIso_cut = 1.8

        muFsrPhotonIdx = [-1] * len(muons)
        eleFsrPhotonIdx = [-1] * len(electrons)
        muPhotondREt2 = [999] * len(muons)
        elePhotondREt2 = [999] * len(electrons)
        fsrPhoton_mydROverEt2 = [-1] * len(fsrPhotons)
        fsrPhoton_myMuonIdx = [-1] * len(fsrPhotons)
        fsrPhoton_myElectronIdx = [-1] * len(fsrPhotons)

        for ifsr, fsr in enumerate(fsrPhotons):
            if fsr.pt < 2 or abs(fsr.eta) > 2.4 or fsr.relIso03 > relIso_cut:
                continue

            dRmin = 0.5
            closestMu=-1
            closestEle=-1

            # Check muons
            for imu, mu in enumerate(muons):
                if not (mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1 and abs(mu.sip3d) < 4 and (mu.isGlobal or (mu.isTracker and mu.nStations > 0))):
                    continue
                dR = deltaR(mu.eta, mu.phi, fsr.eta, fsr.phi)
                if 0.001 < dR < dRmin and dR / (fsr.pt ** 2) < fsr_dRET2_max:
                    dRmin = dR
                    closestMu=imu

            # Check electrons
            for iel, el in enumerate(electrons):
                etaSC = el.eta + el.deltaEtaSC  # still used later, not here
                if not (el.pt > 7 and abs(el.eta) < 2.5 and abs(el.dxy) < 0.5 and abs(el.dz) < 1 and abs(el.sip3d) < 4):
                    continue
                dR = deltaR(el.eta, el.phi, fsr.eta, fsr.phi)
                if 0.001 < dR < dRmin and dR / (fsr.pt ** 2) < fsr_dRET2_max:
                    dRmin = dR
                    closestMu=-1
                    closestEle=iel

            if closestMu >= 0 or closestEle >= 0:
                dREt2 = dRmin / (fsr.pt ** 2)
                fsrPhoton_mydROverEt2[ifsr] = dREt2

                if closestMu >= 0:
                    fsrPhoton_myMuonIdx[ifsr] = closestMu
                    if dREt2 < muPhotondREt2[closestMu]:
                        muPhotondREt2[closestMu] = dREt2
                        muFsrPhotonIdx[closestMu] = ifsr

                if closestEle >= 0:
                    fsrPhoton_myElectronIdx[ifsr] = closestEle
                    if dREt2 < elePhotondREt2[closestEle]:
                        elePhotondREt2[closestEle] = dREt2
                        eleFsrPhotonIdx[closestEle] = ifsr

        # Store results
        event.muFsrPhotonIdx = muFsrPhotonIdx
        event.eleFsrPhotonIdx = eleFsrPhotonIdx

        selectedFSR = []
        for ifsr in muFsrPhotonIdx + eleFsrPhotonIdx :
            if ifsr>=0 : selectedFSR.append(fsrPhotons[ifsr])

        event.allMuons = []
        event.relaxedMuons = []
        event.fullIDMuons = []

        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

    
        for imu, mu in enumerate(muons):
            mu.isRelaxed = False
            mu.isFullID = False
            isoCorr = self.isoFsrCorr(mu, selectedFSR)

            event.allMuons.append(mu)

            if mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1 and abs(mu.sip3d) < 4 and (mu.isGlobal or (mu.isTracker and mu.nStations > 0)):
                mu._wp_Iso = 'LoosePFIso'
                mu.isRelaxed = True
                mu.index = imu
                event.relaxedMuons.append(mu)

                # Full ID with corrected iso
                passMuID = (mu.isPFcand or (mu.highPtId > 0 and mu.pt > 200)) and isoCorr < 0.35
                if passMuID:
                    mu._wp_ID = 'TightID'
                    mu.isFullID = True
                    event.fullIDMuons.append(mu)

        event.allElectrons = []
        event.relaxedElectrons = []
        event.fullIDElectrons = []

        electrons = Collection(event, "Electron")
        fsrPhotons = Collection(event, "FsrPhoton")

        for iel, el in enumerate(electrons):
            el.etaSC = el.eta + el.deltaEtaSC
            el.isRelaxed = False
            el.isFullID = False
            isoCorr = self.isoFsrCorr(el, selectedFSR)
            

            if el.pt > 7 and abs(el.eta) < 2.5 and abs(el.dxy) < 0.5 and abs(el.dz) < 1 and abs(el.sip3d) < 4:
                el._wp_ID = 'wp90iso'
                el.isRelaxed = True
                el.index = iel
                event.relaxedElectrons.append(el)

                # Full ID BDT cuts
                etaSC = abs(el.etaSC)
                mva = el.mvaHZZIso
                if etaSC < 0.8:
                    if el.pt < 10 and mva < 0.9044286167: continue
                    if el.pt >= 10 and mva < 0.1968600840: continue
                elif 0.8 < etaSC < 1.479:
                    if el.pt < 10 and mva < 0.9094166886: continue
                    if el.pt >= 10 and mva < 0.0759172100: continue
                else:
                    if el.pt < 10 and mva < 0.9443653660: continue
                    if el.pt >= 10 and mva < -0.5169136775: continue

                el.isFullID = True
                event.fullIDElectrons.append(el)


    def _select_jets(self, event):

        event.selectedJets = []

        jets = Collection(event, "Jet")
        FsrPhotons = Collection(event, "FsrPhoton")

        for jet in jets:
            if jet.pt <= 20 or abs(jet.eta) >= 2.5:
                continue
            
            jet_isolated = True
            for lep in event.selectedLeptons:
                dR_jet_lep = math.sqrt( (lep.eta - jet.eta)**2 + (lep.phi - jet.phi)**2 ) 
                if dR_jet_lep <= 0.4:
                    jet_isolated = False
            if not jet_isolated:
                continue
                    
            photon_isolated = True
            for photon in FsrPhotons:
                dR_jet_photon = math.sqrt( (photon.eta - jet.eta)**2 + (photon.phi - jet.phi)**2 ) 
                if dR_jet_photon <= 0.4:
                    photon_isolated = False
            if not photon_isolated:
                continue                   
            
            event.selectedJets.append(jet)
        
    def _fill_event_info(self, event):
        out_data = {}
        
        out_data["lumiwgt"] = lumi_dict[self.year]

        ## leptons
        leptons_pt_sorted = sorted(event.selectedLeptons, key=lambda particle: particle.pt, reverse=True)

        lep_pt = []
        lep_eta = []
        lep_phi = []
        lep_pdgId=[]
        el_pt = []
        el_eta = []
        el_phi = []
        el_pdgId =[]
        mu_pt = []
        mu_eta = []
        mu_phi = []
        mu_pdgId = []
        for lep in leptons_pt_sorted:
            lep_pt.append(lep.pt)
            lep_eta.append(lep.eta)
            lep_phi.append(lep.phi)
            lep_pdgId.append(lep.pdgId)
            if abs(lep.pdgId) == 11: # el
                el_pt.append(lep.pt)
                el_eta.append(lep.eta)
                el_phi.append(lep.phi)
                el_pdgId.append(lep.pdgId)
            if abs(lep.pdgId) == 13: # mu
                mu_pt.append(lep.pt)
                mu_eta.append(lep.eta)
                mu_phi.append(lep.phi)
                mu_pdgId.append(lep.pdgId)            
            
        out_data[self.lep_prefix + "pt"] = lep_pt
        out_data[self.lep_prefix + "eta"] = lep_eta
        out_data[self.lep_prefix + "phi"] = lep_phi
        out_data[self.lep_prefix + "pdgId"] = lep_pdgId
        out_data[self.el_prefix + "pt"] = el_pt
        out_data[self.el_prefix + "eta"] = el_eta
        out_data[self.el_prefix + "phi"] = el_phi 
        out_data[self.el_prefix + "pdgId"] = el_pdgId     
        out_data[self.mu_prefix + "pt"] = mu_pt
        out_data[self.mu_prefix + "eta"] = mu_eta
        out_data[self.mu_prefix + "phi"] = mu_phi
        out_data[self.mu_prefix + "pdgId"] = mu_pdgId


        ## Full ID leptons
        fullleptons_pt_sorted = sorted(event.fullIDLeptons, key=lambda particle: particle.pt, reverse=True)

        full_lep_pt = []
        full_lep_eta = []
        full_lep_phi = []
        full_lep_pdgId=[]
        full_el_pt = []
        full_el_eta = []
        full_el_phi = []
        full_el_pdgId =[]
        full_mu_pt = []
        full_mu_eta = []
        full_mu_phi = []
        full_mu_pdgId = []
        for lep in fullleptons_pt_sorted:
            full_lep_pt.append(lep.pt)
            full_lep_eta.append(lep.eta)
            full_lep_phi.append(lep.phi)
            full_lep_pdgId.append(lep.pdgId)
            if abs(lep.pdgId) == 11: # el
                full_el_pt.append(lep.pt)
                full_el_eta.append(lep.eta)
                full_el_phi.append(lep.phi)
                full_el_pdgId.append(lep.pdgId)
            if abs(lep.pdgId) == 13: # mu
                full_mu_pt.append(lep.pt)
                full_mu_eta.append(lep.eta)
                full_mu_phi.append(lep.phi)
                full_mu_pdgId.append(lep.pdgId)            
            
        out_data[self.full_lep_prefix + "pt"] = full_lep_pt
        out_data[self.full_lep_prefix + "eta"] = full_lep_eta
        out_data[self.full_lep_prefix + "phi"] = full_lep_phi
        out_data[self.full_lep_prefix + "pdgId"] = full_lep_pdgId
        out_data[self.full_el_prefix + "pt"] = full_el_pt
        out_data[self.full_el_prefix + "eta"] = full_el_eta
        out_data[self.full_el_prefix + "phi"] = full_el_phi 
        out_data[self.full_el_prefix + "pdgId"] = full_el_pdgId     
        out_data[self.full_mu_prefix + "pt"] = full_mu_pt
        out_data[self.full_mu_prefix + "eta"] = full_mu_eta
        out_data[self.full_mu_prefix + "phi"] = full_mu_phi
        out_data[self.full_mu_prefix + "pdgId"] = full_mu_pdgId
        
        ## jets  
        ak4_pt = []
        ak4_eta = []
        ak4_phi = []
        ak4_mass = []
        ak4_hadronFlavour = []
        
        for jet in event.selectedJets:
            # ak4_bdisc.append(jet.btagDeepFlavB)
            # ak4_cvbdisc.append(jet.btagDeepFlavCvB)
            # ak4_cvldisc.append(jet.btagDeepFlavCvL)
            # ak4_gvudsdisc.append(jet.btagDeepFlavQG)
            ak4_pt.append(jet.pt)
            ak4_eta.append(jet.eta)
            ak4_phi.append(jet.phi)
            ak4_mass.append(jet.mass)
            if self.isMC: ak4_hadronFlavour.append(jet.hadronFlavour)
        
        # out_data[self.jet_prefix + "bdisc"] = ak4_bdisc
        # out_data[self.jet_prefix + "cvbdisc"] = ak4_cvbdisc
        # out_data[self.jet_prefix + "cvldisc"] = ak4_cvldisc 
        # out_data[self.jet_prefix + "gvudsdisc"] = ak4_gvudsdisc 
        out_data[self.jet_prefix + "pt"] = ak4_pt 
        out_data[self.jet_prefix + "eta"] = ak4_eta 
        out_data[self.jet_prefix + "phi"] = ak4_phi 
        out_data[self.jet_prefix + "mass"] = ak4_mass 
        if self.isMC: out_data[self.jet_prefix + "hadronFlavour"] = ak4_hadronFlavour 
           
        ## Z candidates
        Zcandidate_mass = []
        Zcandidate_pt = []
        Zcandidate_eta = []
        Zcandidate_phi = []
        Zcandidate_onshell_mass = []
        Zcandidate_offshell_mass = []
        for Zcandidate in event.ZcandidatesSR:
            Zcandidate_mass.append(Zcandidate.mass)
            Zcandidate_pt.append(Zcandidate.pt)
            Zcandidate_eta.append(Zcandidate.eta)
            Zcandidate_phi.append(Zcandidate.phi)
            if len(event.ZcandidatesSR)>1:
                if Zcandidate.is_onshell: Zcandidate_onshell_mass.append(Zcandidate.mass)
                else: Zcandidate_offshell_mass.append(Zcandidate.mass)
                
        out_data[self.Z_prefix + "mass"] = Zcandidate_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_phi 
        out_data[self.Z_prefix + "onshell_mass"] = Zcandidate_onshell_mass
        out_data[self.Z_prefix + "offshell_mass"] = Zcandidate_offshell_mass

        # ## ZZ candidates
        # # Initialize output lists
        # ZZcandidate_mass = []
        # ZZcandidate_pt = []
        # ZZcandidate_eta = []
        # ZZcandidate_phi = []
        # ZZcandidate_mass_4e=[]
        # ZZcandidate_pt_4e = []
        # ZZcandidate_eta_4e = []
        # ZZcandidate_phi_4e = []
        # ZZcandidate_mass_4mu=[]
        # ZZcandidate_pt_4mu = []
        # ZZcandidate_eta_4mu = []
        # ZZcandidate_phi_4mu = []
        # ZZcandidate_mass_2e2mu=[]
        # ZZcandidate_pt_2e2mu = []
        # ZZcandidate_eta_2e2mu = []
        # ZZcandidate_phi_2e2mu = []

        # for ZZcandidate in event.ZZcandidates:
        #     ZZcandidate_mass.append(ZZcandidate.mass)
        #     ZZcandidate_pt.append(ZZcandidate.pt)
        #     ZZcandidate_eta.append(ZZcandidate.eta)
        #     ZZcandidate_phi.append(ZZcandidate.phi)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZZcandidate.Z1.lep1.pdgId),
        #         abs(ZZcandidate.Z1.lep2.pdgId),
        #         abs(ZZcandidate.Z2.lep1.pdgId),
        #         abs(ZZcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZZcandidate_mass_4e.append(ZZcandidate.mass)
        #         ZZcandidate_pt_4e.append(ZZcandidate.pt)
        #         ZZcandidate_eta_4e.append(ZZcandidate.eta)
        #         ZZcandidate_phi_4e.append(ZZcandidate.phi)
        #     elif lep_ids == {13}:  
        #         ZZcandidate_mass_4mu.append(ZZcandidate.mass)
        #         ZZcandidate_pt_4mu.append(ZZcandidate.pt)
        #         ZZcandidate_eta_4mu.append(ZZcandidate.eta)
        #         ZZcandidate_phi_4mu.append(ZZcandidate.phi)
        #     elif lep_ids == {11, 13}:  
        #         ZZcandidate_mass_2e2mu.append(ZZcandidate.mass)
        #         ZZcandidate_pt_2e2mu.append(ZZcandidate.pt)
        #         ZZcandidate_eta_2e2mu.append(ZZcandidate.eta)
        #         ZZcandidate_phi_2e2mu.append(ZZcandidate.phi)

        # # Store in output dictionary
        # out_data[self.ZZ_prefix + "mass"] = ZZcandidate_mass
        # out_data[self.ZZ_prefix + "pt"] = ZZcandidate_pt
        # out_data[self.ZZ_prefix + "eta"] = ZZcandidate_eta 
        # out_data[self.ZZ_prefix + "phi"] = ZZcandidate_phi

        # out_data[self.ZZ4e_prefix + "mass"] = ZZcandidate_mass_4e
        # out_data[self.ZZ4e_prefix + "pt"] = ZZcandidate_pt_4e
        # out_data[self.ZZ4e_prefix + "eta"] = ZZcandidate_eta_4e
        # out_data[self.ZZ4e_prefix + "phi"] = ZZcandidate_phi_4e

        # out_data[self.ZZ4mu_prefix + "mass"] = ZZcandidate_mass_4mu
        # out_data[self.ZZ4mu_prefix + "pt"] = ZZcandidate_pt_4mu
        # out_data[self.ZZ4mu_prefix + "eta"] = ZZcandidate_eta_4mu
        # out_data[self.ZZ4mu_prefix + "phi"] = ZZcandidate_phi_4mu

        # out_data[self.ZZ2e2mu_prefix + "mass"] =  ZZcandidate_mass_2e2mu
        # out_data[self.ZZ2e2mu_prefix + "pt"] = ZZcandidate_pt_2e2mu
        # out_data[self.ZZ2e2mu_prefix + "eta"] = ZZcandidate_eta_2e2mu
        # out_data[self.ZZ2e2mu_prefix + "phi"] = ZZcandidate_phi_2e2mu

        # # Sanity check 
        # ZZSRcandidate_mass = []
        # ZZSRcandidate_pt = []
        # ZZSRcandidate_eta = []
        # ZZSRcandidate_phi = []
        # ZZSRcandidate_mass_4e=[]
        # ZZSRcandidate_pt_4e = []
        # ZZSRcandidate_eta_4e = []
        # ZZSRcandidate_phi_4e = []
        # ZZSRcandidate_mass_4mu=[]
        # ZZSRcandidate_pt_4mu = []
        # ZZSRcandidate_eta_4mu = []
        # ZZSRcandidate_phi_4mu = []
        # ZZSRcandidate_mass_2e2mu=[]
        # ZZSRcandidate_pt_2e2mu = []
        # ZZSRcandidate_eta_2e2mu = []
        # ZZSRcandidate_phi_2e2mu = []

        # for ZZSRcandidate in event.ZZSRcandidates:
        #     ZZSRcandidate_mass.append(ZZSRcandidate.mass)
        #     ZZSRcandidate_pt.append(ZZSRcandidate.pt)
        #     ZZSRcandidate_eta.append(ZZSRcandidate.eta)
        #     ZZSRcandidate_phi.append(ZZSRcandidate.phi)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZZSRcandidate.Z1.lep1.pdgId),
        #         abs(ZZSRcandidate.Z1.lep2.pdgId),
        #         abs(ZZSRcandidate.Z2.lep1.pdgId),
        #         abs(ZZSRcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZZSRcandidate_mass_4e.append(ZZSRcandidate.mass)
        #         ZZSRcandidate_pt_4e.append(ZZSRcandidate.pt)
        #         ZZSRcandidate_eta_4e.append(ZZSRcandidate.eta)
        #         ZZSRcandidate_phi_4e.append(ZZSRcandidate.phi)
        #     elif lep_ids == {13}:  
        #         ZZSRcandidate_mass_4mu.append(ZZSRcandidate.mass)
        #         ZZSRcandidate_pt_4mu.append(ZZSRcandidate.pt)
        #         ZZSRcandidate_eta_4mu.append(ZZSRcandidate.eta)
        #         ZZSRcandidate_phi_4mu.append(ZZSRcandidate.phi)
        #     elif lep_ids == {11, 13}:  
        #         ZZSRcandidate_mass_2e2mu.append(ZZSRcandidate.mass)
        #         ZZSRcandidate_pt_2e2mu.append(ZZSRcandidate.pt)
        #         ZZSRcandidate_eta_2e2mu.append(ZZSRcandidate.eta)
        #         ZZSRcandidate_phi_2e2mu.append(ZZSRcandidate.phi)

        # # Store in output dictionary
        # out_data[self.ZZSR_prefix + "mass"] = ZZSRcandidate_mass
        # out_data[self.ZZSR_prefix + "pt"] = ZZSRcandidate_pt
        # out_data[self.ZZSR_prefix + "eta"] = ZZSRcandidate_eta 
        # out_data[self.ZZSR_prefix + "phi"] = ZZSRcandidate_phi

        # out_data[self.ZZSR4e_prefix + "mass"] = ZZSRcandidate_mass_4e
        # out_data[self.ZZSR4e_prefix + "pt"] = ZZSRcandidate_pt_4e
        # out_data[self.ZZSR4e_prefix + "eta"] = ZZSRcandidate_eta_4e
        # out_data[self.ZZSR4e_prefix + "phi"] = ZZSRcandidate_phi_4e

        # out_data[self.ZZSR4mu_prefix + "mass"] = ZZSRcandidate_mass_4mu
        # out_data[self.ZZSR4mu_prefix + "pt"] = ZZSRcandidate_pt_4mu
        # out_data[self.ZZSR4mu_prefix + "eta"] = ZZSRcandidate_eta_4mu
        # out_data[self.ZZSR4mu_prefix + "phi"] = ZZSRcandidate_phi_4mu

        # out_data[self.ZZSR2e2mu_prefix + "mass"] =  ZZSRcandidate_mass_2e2mu
        # out_data[self.ZZSR2e2mu_prefix + "pt"] = ZZSRcandidate_pt_2e2mu
        # out_data[self.ZZSR2e2mu_prefix + "eta"] = ZZSRcandidate_eta_2e2mu
        # out_data[self.ZZSR2e2mu_prefix + "phi"] = ZZSRcandidate_phi_2e2mu

        # # Similar structure for Higgs candidates
        # Hcandidate_mass = []
        # Hcandidate_pt = []
        # Hcandidate_eta = []
        # Hcandidate_phi = []
        # Hcandidate_mass_4e=[]
        # Hcandidate_pt_4e = []
        # Hcandidate_eta_4e = []
        # Hcandidate_phi_4e = []
        # Hcandidate_mass_4mu=[]
        # Hcandidate_pt_4mu = []
        # Hcandidate_eta_4mu = []
        # Hcandidate_phi_4mu = []
        # Hcandidate_mass_2e2mu=[]
        # Hcandidate_pt_2e2mu = []
        # Hcandidate_eta_2e2mu = []
        # Hcandidate_phi_2e2mu = []

        # for Hcandidate in event.Hcandidates:
        #     Hcandidate_mass.append(Hcandidate.mass)
        #     Hcandidate_pt.append(Hcandidate.pt)
        #     Hcandidate_eta.append(Hcandidate.eta)
        #     Hcandidate_phi.append(Hcandidate.phi)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(Hcandidate.Z1.lep1.pdgId),
        #         abs(Hcandidate.Z1.lep2.pdgId),
        #         abs(Hcandidate.Z2.lep1.pdgId),
        #         abs(Hcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         Hcandidate_mass_4e.append(Hcandidate.mass)
        #         Hcandidate_pt_4e.append(Hcandidate.pt)
        #         Hcandidate_eta_4e.append(Hcandidate.eta)
        #         Hcandidate_phi_4e.append(Hcandidate.phi)
        #     elif lep_ids == {13}:  
        #         Hcandidate_mass_4mu.append(Hcandidate.mass)
        #         Hcandidate_pt_4mu.append(Hcandidate.pt)
        #         Hcandidate_eta_4mu.append(Hcandidate.eta)
        #         Hcandidate_phi_4mu.append(Hcandidate.phi)
        #     elif lep_ids == {11, 13}:  
        #         Hcandidate_mass_2e2mu.append(Hcandidate.mass)
        #         Hcandidate_pt_2e2mu.append(Hcandidate.pt)
        #         Hcandidate_eta_2e2mu.append(Hcandidate.eta)
        #         Hcandidate_phi_2e2mu.append(Hcandidate.phi)
                
        # out_data[self.H_prefix + "mass"] = Hcandidate_mass
        # out_data[self.H_prefix + "pt"] = Hcandidate_pt
        # out_data[self.H_prefix + "eta"] = Hcandidate_eta 
        # out_data[self.H_prefix + "phi"] = Hcandidate_phi

        # out_data[self.H4e_prefix + "mass"] = Hcandidate_mass_4e
        # out_data[self.H4e_prefix + "pt"] = Hcandidate_pt_4e
        # out_data[self.H4e_prefix + "eta"] = Hcandidate_eta_4e
        # out_data[self.H4e_prefix + "phi"] = Hcandidate_phi_4e

        # out_data[self.H4mu_prefix + "mass"] = Hcandidate_mass_4mu
        # out_data[self.H4mu_prefix + "pt"] = Hcandidate_pt_4mu
        # out_data[self.H4mu_prefix + "eta"] = Hcandidate_eta_4mu
        # out_data[self.H4mu_prefix + "phi"] = Hcandidate_phi_4mu

        # out_data[self.H2e2mu_prefix + "mass"] =  Hcandidate_mass_2e2mu
        # out_data[self.H2e2mu_prefix + "pt"] = Hcandidate_pt_2e2mu
        # out_data[self.H2e2mu_prefix + "eta"] = Hcandidate_eta_2e2mu
        # out_data[self.H2e2mu_prefix + "phi"] = Hcandidate_phi_2e2mu
                
        # ## ZLL candidates
        # # Initialize output lists
        # ZLLcandidate_mass = []
        # ZLLcandidate_pt = []
        # ZLLcandidate_eta = []
        # ZLLcandidate_phi = []
        # ZLLcandidate_mass_4e=[]
        # ZLLcandidate_pt_4e = []
        # ZLLcandidate_eta_4e = []
        # ZLLcandidate_phi_4e = []
        # ZLLcandidate_lep3_pt_4e = []
        # ZLLcandidate_lep3_eta_4e = []
        # ZLLcandidate_lep3_pdgId_4e = []
        # ZLLcandidate_lep4_pt_4e = []
        # ZLLcandidate_lep4_eta_4e = []
        # ZLLcandidate_lep4_pdgId_4e = []
        # ZLLcandidate_lep3_lostHits_4e = []
        # ZLLcandidate_lep4_lostHits_4e = []
        # ZLLcandidate_lep3_isGood_4e = []
        # ZLLcandidate_lep4_isGood_4e = []
        # ZLLcandidate_mass_4mu=[]
        # ZLLcandidate_pt_4mu = []
        # ZLLcandidate_eta_4mu = []
        # ZLLcandidate_phi_4mu = []
        # ZLLcandidate_mass_2e2mu=[]
        # ZLLcandidate_pt_2e2mu = []
        # ZLLcandidate_eta_2e2mu = []
        # ZLLcandidate_phi_2e2mu = []

        # for ZLLcandidate in event.ZLLcandidates:
        #     ZLLcandidate_mass.append(ZLLcandidate.mass)
        #     ZLLcandidate_pt.append(ZLLcandidate.pt)
        #     ZLLcandidate_eta.append(ZLLcandidate.eta)
        #     ZLLcandidate_phi.append(ZLLcandidate.phi)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZLLcandidate.Z1.lep1.pdgId),
        #         abs(ZLLcandidate.Z1.lep2.pdgId),
        #         abs(ZLLcandidate.Z2.lep1.pdgId),
        #         abs(ZLLcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZLLcandidate_mass_4e.append(ZLLcandidate.mass)
        #         ZLLcandidate_pt_4e.append(ZLLcandidate.pt)
        #         ZLLcandidate_eta_4e.append(ZLLcandidate.eta)
        #         ZLLcandidate_phi_4e.append(ZLLcandidate.phi)
        #         ZLLcandidate_lep3_pt_4e.append(ZLLcandidate.Z2.lep1.pt)
        #         ZLLcandidate_lep3_eta_4e.append(ZLLcandidate.Z2.lep1.eta)
        #         ZLLcandidate_lep3_pdgId_4e.append(ZLLcandidate.Z2.lep1.pdgId)
        #         ZLLcandidate_lep4_pt_4e.append(ZLLcandidate.Z2.lep2.pt)
        #         ZLLcandidate_lep4_eta_4e.append(ZLLcandidate.Z2.lep2.eta)
        #         ZLLcandidate_lep4_pdgId_4e.append(ZLLcandidate.Z2.lep2.pdgId)
        #         ZLLcandidate_lep3_lostHits_4e.append(ZLLcandidate.Z2.lep1.lostHits)
        #         ZLLcandidate_lep4_lostHits_4e.append(ZLLcandidate.Z2.lep2.lostHits)
        #         ZLLcandidate_lep3_isGood_4e.append(ZLLcandidate.Z2.lep1.isFullID)
        #         ZLLcandidate_lep4_isGood_4e.append(ZLLcandidate.Z2.lep2.isFullID)
        #     elif lep_ids == {13}:  
        #         ZLLcandidate_mass_4mu.append(ZLLcandidate.mass)
        #         ZLLcandidate_pt_4mu.append(ZLLcandidate.pt)
        #         ZLLcandidate_eta_4mu.append(ZLLcandidate.eta)
        #         ZLLcandidate_phi_4mu.append(ZLLcandidate.phi)
        #     elif lep_ids == {11, 13}:  
        #         ZLLcandidate_mass_2e2mu.append(ZLLcandidate.mass)
        #         ZLLcandidate_pt_2e2mu.append(ZLLcandidate.pt)
        #         ZLLcandidate_eta_2e2mu.append(ZLLcandidate.eta)
        #         ZLLcandidate_phi_2e2mu.append(ZLLcandidate.phi)

        # # Store in output dictionary
        # out_data[self.ZLL_prefix + "mass"] = ZLLcandidate_mass
        # out_data[self.ZLL_prefix + "pt"] = ZLLcandidate_pt
        # out_data[self.ZLL_prefix + "eta"] = ZLLcandidate_eta 
        # out_data[self.ZLL_prefix + "phi"] = ZLLcandidate_phi

        # out_data[self.ZLL4e_prefix + "mass"] = ZLLcandidate_mass_4e
        # out_data[self.ZLL4e_prefix + "pt"] = ZLLcandidate_pt_4e
        # out_data[self.ZLL4e_prefix + "eta"] = ZLLcandidate_eta_4e
        # out_data[self.ZLL4e_prefix + "phi"] = ZLLcandidate_phi_4e
        # out_data[self.ZLL4e_prefix + "lep3_pt"] = ZLLcandidate_lep3_pt_4e
        # out_data[self.ZLL4e_prefix + "lep3_eta"] = ZLLcandidate_lep3_eta_4e
        # out_data[self.ZLL4e_prefix + "lep3_pdgId"] = ZLLcandidate_lep3_pdgId_4e
        # out_data[self.ZLL4e_prefix + "lep4_pt"] = ZLLcandidate_lep4_pt_4e
        # out_data[self.ZLL4e_prefix + "lep4_eta"] = ZLLcandidate_lep4_eta_4e
        # out_data[self.ZLL4e_prefix + "lep4_pdgId"] = ZLLcandidate_lep4_pdgId_4e
        # out_data[self.ZLL4e_prefix + "lep3_lostHits"] = ZLLcandidate_lep3_lostHits_4e
        # out_data[self.ZLL4e_prefix + "lep4_lostHits"] = ZLLcandidate_lep4_lostHits_4e
        # out_data[self.ZLL4e_prefix + "lep3_isGood"] = ZLLcandidate_lep3_isGood_4e
        # out_data[self.ZLL4e_prefix + "lep4_isGood"] = ZLLcandidate_lep4_isGood_4e

        # out_data[self.ZLL4mu_prefix + "mass"] = ZLLcandidate_mass_4mu
        # out_data[self.ZLL4mu_prefix + "pt"] = ZLLcandidate_pt_4mu
        # out_data[self.ZLL4mu_prefix + "eta"] = ZLLcandidate_eta_4mu
        # out_data[self.ZLL4mu_prefix + "phi"] = ZLLcandidate_phi_4mu

        # out_data[self.ZLL2e2mu_prefix + "mass"] =  ZLLcandidate_mass_2e2mu
        # out_data[self.ZLL2e2mu_prefix + "pt"] = ZLLcandidate_pt_2e2mu
        # out_data[self.ZLL2e2mu_prefix + "eta"] = ZLLcandidate_eta_2e2mu
        # out_data[self.ZLL2e2mu_prefix + "phi"] = ZLLcandidate_phi_2e2mu

        # # 2P2F
        # ZLL2P2Fcandidate_mass = []
        # ZLL2P2Fcandidate_pt = []
        # ZLL2P2Fcandidate_eta = []
        # ZLL2P2Fcandidate_phi = []
        # ZLL2P2Fcandidate_lep3_pt = []
        # ZLL2P2Fcandidate_lep3_eta = []
        # ZLL2P2Fcandidate_lep3_pdgId = []
        # ZLL2P2Fcandidate_lep4_pt = []
        # ZLL2P2Fcandidate_lep4_eta = []
        # ZLL2P2Fcandidate_lep4_pdgId = []
        # ZLL2P2Fcandidate_mass_4e=[]
        # ZLL2P2Fcandidate_pt_4e = []
        # ZLL2P2Fcandidate_eta_4e = []
        # ZLL2P2Fcandidate_phi_4e = []
        # ZLL2P2Fcandidate_lep3_pt_4e = []
        # ZLL2P2Fcandidate_lep3_eta_4e = []
        # ZLL2P2Fcandidate_lep3_pdgId_4e = []
        # ZLL2P2Fcandidate_lep4_pt_4e = []
        # ZLL2P2Fcandidate_lep4_eta_4e = []
        # ZLL2P2Fcandidate_lep4_pdgId_4e = []
        # ZLL2P2Fcandidate_mass_4mu=[]
        # ZLL2P2Fcandidate_pt_4mu = []
        # ZLL2P2Fcandidate_eta_4mu = []
        # ZLL2P2Fcandidate_phi_4mu = []
        # ZLL2P2Fcandidate_lep3_pt_4mu = []
        # ZLL2P2Fcandidate_lep3_eta_4mu = []
        # ZLL2P2Fcandidate_lep3_pdgId_4mu = []
        # ZLL2P2Fcandidate_lep4_pt_4mu = []
        # ZLL2P2Fcandidate_lep4_eta_4mu = []
        # ZLL2P2Fcandidate_lep4_pdgId_4mu = []
        # ZLL2P2Fcandidate_mass_2e2mu=[]
        # ZLL2P2Fcandidate_pt_2e2mu = []
        # ZLL2P2Fcandidate_eta_2e2mu = []
        # ZLL2P2Fcandidate_phi_2e2mu = []
        # ZLL2P2Fcandidate_lep3_pt_2e2mu = []
        # ZLL2P2Fcandidate_lep3_eta_2e2mu = []
        # ZLL2P2Fcandidate_lep3_pdgId_2e2mu = []
        # ZLL2P2Fcandidate_lep4_pt_2e2mu = []
        # ZLL2P2Fcandidate_lep4_eta_2e2mu = []
        # ZLL2P2Fcandidate_lep4_pdgId_2e2mu = []

        # for ZLL2P2Fcandidate in event.ZLL2P2Fcandidates:
        #     ZLL2P2Fcandidate_mass.append(ZLL2P2Fcandidate.mass)
        #     ZLL2P2Fcandidate_pt.append(ZLL2P2Fcandidate.pt)
        #     ZLL2P2Fcandidate_eta.append(ZLL2P2Fcandidate.eta)
        #     ZLL2P2Fcandidate_phi.append(ZLL2P2Fcandidate.phi)
        #     ZLL2P2Fcandidate_lep3_pt.append(ZLL2P2Fcandidate.Z2.lep1.pt)
        #     ZLL2P2Fcandidate_lep3_eta.append(ZLL2P2Fcandidate.Z2.lep1.eta)
        #     ZLL2P2Fcandidate_lep3_pdgId.append(ZLL2P2Fcandidate.Z2.lep1.pdgId)
        #     ZLL2P2Fcandidate_lep4_pt.append(ZLL2P2Fcandidate.Z2.lep2.pt)
        #     ZLL2P2Fcandidate_lep4_eta.append(ZLL2P2Fcandidate.Z2.lep2.eta)
        #     ZLL2P2Fcandidate_lep4_pdgId.append(ZLL2P2Fcandidate.Z2.lep2.pdgId)
        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZLL2P2Fcandidate.Z1.lep1.pdgId),
        #         abs(ZLL2P2Fcandidate.Z1.lep2.pdgId),
        #         abs(ZLL2P2Fcandidate.Z2.lep1.pdgId),
        #         abs(ZLL2P2Fcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZLL2P2Fcandidate_mass_4e.append(ZLL2P2Fcandidate.mass)
        #         ZLL2P2Fcandidate_pt_4e.append(ZLL2P2Fcandidate.pt)
        #         ZLL2P2Fcandidate_eta_4e.append(ZLL2P2Fcandidate.eta)
        #         ZLL2P2Fcandidate_phi_4e.append(ZLL2P2Fcandidate.phi)
        #         ZLL2P2Fcandidate_lep3_pt_4e.append(ZLL2P2Fcandidate.Z2.lep1.pt)
        #         ZLL2P2Fcandidate_lep3_eta_4e.append(ZLL2P2Fcandidate.Z2.lep1.eta)
        #         ZLL2P2Fcandidate_lep3_pdgId_4e.append(ZLL2P2Fcandidate.Z2.lep1.pdgId)
        #         ZLL2P2Fcandidate_lep4_pt_4e.append(ZLL2P2Fcandidate.Z2.lep2.pt)
        #         ZLL2P2Fcandidate_lep4_eta_4e.append(ZLL2P2Fcandidate.Z2.lep2.eta)
        #         ZLL2P2Fcandidate_lep4_pdgId_4e.append(ZLL2P2Fcandidate.Z2.lep2.pdgId)
        #     elif lep_ids == {13}:  
        #         ZLL2P2Fcandidate_mass_4mu.append(ZLL2P2Fcandidate.mass)
        #         ZLL2P2Fcandidate_pt_4mu.append(ZLL2P2Fcandidate.pt)
        #         ZLL2P2Fcandidate_eta_4mu.append(ZLL2P2Fcandidate.eta)
        #         ZLL2P2Fcandidate_phi_4mu.append(ZLL2P2Fcandidate.phi)
        #         ZLL2P2Fcandidate_lep3_pt_4mu.append(ZLL2P2Fcandidate.Z2.lep1.pt)
        #         ZLL2P2Fcandidate_lep3_eta_4mu.append(ZLL2P2Fcandidate.Z2.lep1.eta)
        #         ZLL2P2Fcandidate_lep3_pdgId_4mu.append(ZLL2P2Fcandidate.Z2.lep1.pdgId)
        #         ZLL2P2Fcandidate_lep4_pt_4mu.append(ZLL2P2Fcandidate.Z2.lep2.pt)
        #         ZLL2P2Fcandidate_lep4_eta_4mu.append(ZLL2P2Fcandidate.Z2.lep2.eta)
        #         ZLL2P2Fcandidate_lep4_pdgId_4mu.append(ZLL2P2Fcandidate.Z2.lep2.pdgId)
        #     elif lep_ids == {11, 13}:  
        #         ZLL2P2Fcandidate_mass_2e2mu.append(ZLL2P2Fcandidate.mass)
        #         ZLL2P2Fcandidate_pt_2e2mu.append(ZLL2P2Fcandidate.pt)
        #         ZLL2P2Fcandidate_eta_2e2mu.append(ZLL2P2Fcandidate.eta)
        #         ZLL2P2Fcandidate_phi_2e2mu.append(ZLL2P2Fcandidate.phi)
        #         ZLL2P2Fcandidate_lep3_pt_2e2mu.append(ZLL2P2Fcandidate.Z2.lep1.pt)
        #         ZLL2P2Fcandidate_lep3_eta_2e2mu.append(ZLL2P2Fcandidate.Z2.lep1.eta)
        #         ZLL2P2Fcandidate_lep3_pdgId_2e2mu.append(ZLL2P2Fcandidate.Z2.lep1.pdgId)
        #         ZLL2P2Fcandidate_lep4_pt_2e2mu.append(ZLL2P2Fcandidate.Z2.lep2.pt)
        #         ZLL2P2Fcandidate_lep4_eta_2e2mu.append(ZLL2P2Fcandidate.Z2.lep2.eta)
        #         ZLL2P2Fcandidate_lep4_pdgId_2e2mu.append(ZLL2P2Fcandidate.Z2.lep2.pdgId)

        # # Store in output dictionary
        # out_data[self.ZLL2P2F_prefix + "mass"] = ZLL2P2Fcandidate_mass
        # out_data[self.ZLL2P2F_prefix + "pt"] = ZLL2P2Fcandidate_pt
        # out_data[self.ZLL2P2F_prefix + "eta"] = ZLL2P2Fcandidate_eta 
        # out_data[self.ZLL2P2F_prefix + "phi"] = ZLL2P2Fcandidate_phi
        # out_data[self.ZLL2P2F_prefix + "lep3_pt"] = ZLL2P2Fcandidate_lep3_pt
        # out_data[self.ZLL2P2F_prefix + "lep3_eta"] = ZLL2P2Fcandidate_lep3_eta
        # out_data[self.ZLL2P2F_prefix + "lep3_pdgId"] = ZLL2P2Fcandidate_lep3_pdgId
        # out_data[self.ZLL2P2F_prefix + "lep4_pt"] = ZLL2P2Fcandidate_lep4_pt
        # out_data[self.ZLL2P2F_prefix + "lep4_eta"] = ZLL2P2Fcandidate_lep4_eta
        # out_data[self.ZLL2P2F_prefix + "lep4_pdgId"] = ZLL2P2Fcandidate_lep4_pdgId


        # out_data[self.ZLL2P2F4e_prefix + "mass"] = ZLL2P2Fcandidate_mass_4e
        # out_data[self.ZLL2P2F4e_prefix + "pt"] = ZLL2P2Fcandidate_pt_4e
        # out_data[self.ZLL2P2F4e_prefix + "eta"] = ZLL2P2Fcandidate_eta_4e
        # out_data[self.ZLL2P2F4e_prefix + "phi"] = ZLL2P2Fcandidate_phi_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep3_pt"] = ZLL2P2Fcandidate_lep3_pt_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep3_eta"] = ZLL2P2Fcandidate_lep3_eta_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep3_pdgId"] = ZLL2P2Fcandidate_lep3_pdgId_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep4_pt"] = ZLL2P2Fcandidate_lep4_pt_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep4_eta"] = ZLL2P2Fcandidate_lep4_eta_4e
        # out_data[self.ZLL2P2F4e_prefix + "lep4_pdgId"] = ZLL2P2Fcandidate_lep4_pdgId_4e

        # out_data[self.ZLL2P2F4mu_prefix + "mass"] = ZLL2P2Fcandidate_mass_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "pt"] = ZLL2P2Fcandidate_pt_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "eta"] = ZLL2P2Fcandidate_eta_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "phi"] = ZLL2P2Fcandidate_phi_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep3_pt"] = ZLL2P2Fcandidate_lep3_pt_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep3_eta"] = ZLL2P2Fcandidate_lep3_eta_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep3_pdgId"] = ZLL2P2Fcandidate_lep3_pdgId_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep4_pt"] = ZLL2P2Fcandidate_lep4_pt_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep4_eta"] = ZLL2P2Fcandidate_lep4_eta_4mu
        # out_data[self.ZLL2P2F4mu_prefix + "lep4_pdgId"] = ZLL2P2Fcandidate_lep4_pdgId_4mu

        # out_data[self.ZLL2P2F2e2mu_prefix + "mass"] =  ZLL2P2Fcandidate_mass_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "pt"] = ZLL2P2Fcandidate_pt_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "eta"] = ZLL2P2Fcandidate_eta_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "phi"] = ZLL2P2Fcandidate_phi_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep3_pt"] = ZLL2P2Fcandidate_lep3_pt_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep3_eta"] = ZLL2P2Fcandidate_lep3_eta_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep3_pdgId"] = ZLL2P2Fcandidate_lep3_pdgId_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep4_pt"] = ZLL2P2Fcandidate_lep4_pt_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep4_eta"] = ZLL2P2Fcandidate_lep4_eta_2e2mu
        # out_data[self.ZLL2P2F2e2mu_prefix + "lep4_pdgId"] = ZLL2P2Fcandidate_lep4_pdgId_2e2mu

        # # 3P1F
        # ZLL3P1Fcandidate_mass = []
        # ZLL3P1Fcandidate_pt = []
        # ZLL3P1Fcandidate_eta = []
        # ZLL3P1Fcandidate_phi = []
        # ZLL3P1Fcandidate_lep3_pt = []
        # ZLL3P1Fcandidate_lep3_eta = []
        # ZLL3P1Fcandidate_lep3_pdgId = []
        # ZLL3P1Fcandidate_lep4_pt = []
        # ZLL3P1Fcandidate_lep4_eta = []
        # ZLL3P1Fcandidate_lep4_pdgId = []
        # ZLL3P1Fcandidate_mass_4e=[]
        # ZLL3P1Fcandidate_pt_4e = []
        # ZLL3P1Fcandidate_eta_4e = []
        # ZLL3P1Fcandidate_phi_4e = []
        # ZLL3P1Fcandidate_lep3_pt_4e = []
        # ZLL3P1Fcandidate_lep3_eta_4e = []
        # ZLL3P1Fcandidate_lep3_pdgId_4e = []
        # ZLL3P1Fcandidate_lep4_pt_4e = []
        # ZLL3P1Fcandidate_lep4_eta_4e = []
        # ZLL3P1Fcandidate_lep4_pdgId_4e = []
        # ZLL3P1Fcandidate_mass_4mu=[]
        # ZLL3P1Fcandidate_pt_4mu = []
        # ZLL3P1Fcandidate_eta_4mu = []
        # ZLL3P1Fcandidate_phi_4mu = []
        # ZLL3P1Fcandidate_lep3_pt_4mu = []
        # ZLL3P1Fcandidate_lep3_eta_4mu = []
        # ZLL3P1Fcandidate_lep3_pdgId_4mu = []
        # ZLL3P1Fcandidate_lep4_pt_4mu = []
        # ZLL3P1Fcandidate_lep4_eta_4mu = []
        # ZLL3P1Fcandidate_lep4_pdgId_4mu = []
        # ZLL3P1Fcandidate_mass_2e2mu=[]
        # ZLL3P1Fcandidate_pt_2e2mu = []
        # ZLL3P1Fcandidate_eta_2e2mu = []
        # ZLL3P1Fcandidate_phi_2e2mu = []
        # ZLL3P1Fcandidate_lep3_pt_2e2mu = []
        # ZLL3P1Fcandidate_lep3_eta_2e2mu = []
        # ZLL3P1Fcandidate_lep3_pdgId_2e2mu = []
        # ZLL3P1Fcandidate_lep4_pt_2e2mu = []
        # ZLL3P1Fcandidate_lep4_eta_2e2mu = []
        # ZLL3P1Fcandidate_lep4_pdgId_2e2mu = []

        # for ZLL3P1Fcandidate in event.ZLL3P1Fcandidates:
        #     ZLL3P1Fcandidate_mass.append(ZLL3P1Fcandidate.mass)
        #     ZLL3P1Fcandidate_pt.append(ZLL3P1Fcandidate.pt)
        #     ZLL3P1Fcandidate_eta.append(ZLL3P1Fcandidate.eta)
        #     ZLL3P1Fcandidate_phi.append(ZLL3P1Fcandidate.phi)
        #     ZLL3P1Fcandidate_lep3_pt.append(ZLL3P1Fcandidate.Z2.lep1.pt)
        #     ZLL3P1Fcandidate_lep3_eta.append(ZLL3P1Fcandidate.Z2.lep1.eta)
        #     ZLL3P1Fcandidate_lep3_pdgId.append(ZLL3P1Fcandidate.Z2.lep1.pdgId)
        #     ZLL3P1Fcandidate_lep4_pt.append(ZLL3P1Fcandidate.Z2.lep2.pt)
        #     ZLL3P1Fcandidate_lep4_eta.append(ZLL3P1Fcandidate.Z2.lep2.eta)
        #     ZLL3P1Fcandidate_lep4_pdgId.append(ZLL3P1Fcandidate.Z2.lep2.pdgId)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZLL3P1Fcandidate.Z1.lep1.pdgId),
        #         abs(ZLL3P1Fcandidate.Z1.lep2.pdgId),
        #         abs(ZLL3P1Fcandidate.Z2.lep1.pdgId),
        #         abs(ZLL3P1Fcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZLL3P1Fcandidate_mass_4e.append(ZLL3P1Fcandidate.mass)
        #         ZLL3P1Fcandidate_pt_4e.append(ZLL3P1Fcandidate.pt)
        #         ZLL3P1Fcandidate_eta_4e.append(ZLL3P1Fcandidate.eta)
        #         ZLL3P1Fcandidate_phi_4e.append(ZLL3P1Fcandidate.phi)
        #         ZLL3P1Fcandidate_lep3_pt_4e.append(ZLL3P1Fcandidate.Z2.lep1.pt)
        #         ZLL3P1Fcandidate_lep3_eta_4e.append(ZLL3P1Fcandidate.Z2.lep1.eta)
        #         ZLL3P1Fcandidate_lep3_pdgId_4e.append(ZLL3P1Fcandidate.Z2.lep1.pdgId)
        #         ZLL3P1Fcandidate_lep4_pt_4e.append(ZLL3P1Fcandidate.Z2.lep2.pt)
        #         ZLL3P1Fcandidate_lep4_eta_4e.append(ZLL3P1Fcandidate.Z2.lep2.eta)
        #         ZLL3P1Fcandidate_lep4_pdgId_4e.append(ZLL3P1Fcandidate.Z2.lep2.pdgId)
        #     elif lep_ids == {13}:  
        #         ZLL3P1Fcandidate_mass_4mu.append(ZLL3P1Fcandidate.mass)
        #         ZLL3P1Fcandidate_pt_4mu.append(ZLL3P1Fcandidate.pt)
        #         ZLL3P1Fcandidate_eta_4mu.append(ZLL3P1Fcandidate.eta)
        #         ZLL3P1Fcandidate_phi_4mu.append(ZLL3P1Fcandidate.phi)
        #         ZLL3P1Fcandidate_lep3_pt_4mu.append(ZLL3P1Fcandidate.Z2.lep1.pt)
        #         ZLL3P1Fcandidate_lep3_eta_4mu.append(ZLL3P1Fcandidate.Z2.lep1.eta)
        #         ZLL3P1Fcandidate_lep3_pdgId_4mu.append(ZLL3P1Fcandidate.Z2.lep1.pdgId)
        #         ZLL3P1Fcandidate_lep4_pt_4mu.append(ZLL3P1Fcandidate.Z2.lep2.pt)
        #         ZLL3P1Fcandidate_lep4_eta_4mu.append(ZLL3P1Fcandidate.Z2.lep2.eta)
        #         ZLL3P1Fcandidate_lep4_pdgId_4mu.append(ZLL3P1Fcandidate.Z2.lep2.pdgId)
        #     elif lep_ids == {11, 13}:  
        #         ZLL3P1Fcandidate_mass_2e2mu.append(ZLL3P1Fcandidate.mass)
        #         ZLL3P1Fcandidate_pt_2e2mu.append(ZLL3P1Fcandidate.pt)
        #         ZLL3P1Fcandidate_eta_2e2mu.append(ZLL3P1Fcandidate.eta)
        #         ZLL3P1Fcandidate_phi_2e2mu.append(ZLL3P1Fcandidate.phi)
        #         ZLL3P1Fcandidate_lep3_pt_2e2mu.append(ZLL3P1Fcandidate.Z2.lep1.pt)
        #         ZLL3P1Fcandidate_lep3_eta_2e2mu.append(ZLL3P1Fcandidate.Z2.lep1.eta)
        #         ZLL3P1Fcandidate_lep3_pdgId_2e2mu.append(ZLL3P1Fcandidate.Z2.lep1.pdgId)
        #         ZLL3P1Fcandidate_lep4_pt_2e2mu.append(ZLL3P1Fcandidate.Z2.lep2.pt)
        #         ZLL3P1Fcandidate_lep4_eta_2e2mu.append(ZLL3P1Fcandidate.Z2.lep2.eta)
        #         ZLL3P1Fcandidate_lep4_pdgId_2e2mu.append(ZLL3P1Fcandidate.Z2.lep2.pdgId)

        # # Store in output dictionary
        # out_data[self.ZLL3P1F_prefix + "mass"] = ZLL3P1Fcandidate_mass
        # out_data[self.ZLL3P1F_prefix + "pt"] = ZLL3P1Fcandidate_pt
        # out_data[self.ZLL3P1F_prefix + "eta"] = ZLL3P1Fcandidate_eta 
        # out_data[self.ZLL3P1F_prefix + "phi"] = ZLL3P1Fcandidate_phi
        # out_data[self.ZLL3P1F_prefix + "lep3_pt"] = ZLL3P1Fcandidate_lep3_pt
        # out_data[self.ZLL3P1F_prefix + "lep3_eta"] = ZLL3P1Fcandidate_lep3_eta
        # out_data[self.ZLL3P1F_prefix + "lep3_pdgId"] = ZLL3P1Fcandidate_lep3_pdgId
        # out_data[self.ZLL3P1F_prefix + "lep4_pt"] = ZLL3P1Fcandidate_lep4_pt
        # out_data[self.ZLL3P1F_prefix + "lep4_eta"] = ZLL3P1Fcandidate_lep4_eta
        # out_data[self.ZLL3P1F_prefix + "lep4_pdgId"] = ZLL3P1Fcandidate_lep4_pdgId

        # out_data[self.ZLL3P1F4e_prefix + "mass"] = ZLL3P1Fcandidate_mass_4e
        # out_data[self.ZLL3P1F4e_prefix + "pt"] = ZLL3P1Fcandidate_pt_4e
        # out_data[self.ZLL3P1F4e_prefix + "eta"] = ZLL3P1Fcandidate_eta_4e
        # out_data[self.ZLL3P1F4e_prefix + "phi"] = ZLL3P1Fcandidate_phi_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep3_pt"] = ZLL3P1Fcandidate_lep3_pt_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep3_eta"] = ZLL3P1Fcandidate_lep3_eta_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep3_pdgId"] = ZLL3P1Fcandidate_lep3_pdgId_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep4_pt"] = ZLL3P1Fcandidate_lep4_pt_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep4_eta"] = ZLL3P1Fcandidate_lep4_eta_4e
        # out_data[self.ZLL3P1F4e_prefix + "lep4_pdgId"] = ZLL3P1Fcandidate_lep4_pdgId_4e

        # out_data[self.ZLL3P1F4mu_prefix + "mass"] = ZLL3P1Fcandidate_mass_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "pt"] = ZLL3P1Fcandidate_pt_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "eta"] = ZLL3P1Fcandidate_eta_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "phi"] = ZLL3P1Fcandidate_phi_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep3_pt"] = ZLL3P1Fcandidate_lep3_pt_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep3_eta"] = ZLL3P1Fcandidate_lep3_eta_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep3_pdgId"] = ZLL3P1Fcandidate_lep3_pdgId_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep4_pt"] = ZLL3P1Fcandidate_lep4_pt_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep4_eta"] = ZLL3P1Fcandidate_lep4_eta_4mu
        # out_data[self.ZLL3P1F4mu_prefix + "lep4_pdgId"] = ZLL3P1Fcandidate_lep4_pdgId_4mu

        # out_data[self.ZLL3P1F2e2mu_prefix + "mass"] =  ZLL3P1Fcandidate_mass_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "pt"] = ZLL3P1Fcandidate_pt_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "eta"] = ZLL3P1Fcandidate_eta_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "phi"] = ZLL3P1Fcandidate_phi_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep3_pt"] = ZLL3P1Fcandidate_lep3_pt_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep3_eta"] = ZLL3P1Fcandidate_lep3_eta_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep3_pdgId"] = ZLL3P1Fcandidate_lep3_pdgId_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep4_pt"] = ZLL3P1Fcandidate_lep4_pt_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep4_eta"] = ZLL3P1Fcandidate_lep4_eta_2e2mu
        # out_data[self.ZLL3P1F2e2mu_prefix + "lep4_pdgId"] = ZLL3P1Fcandidate_lep4_pdgId_2e2mu

        # # SSCR
        # ZLLSSCRcandidate_mass = []
        # ZLLSSCRcandidate_pt = []
        # ZLLSSCRcandidate_eta = []
        # ZLLSSCRcandidate_phi = []

        # ZLLSSCRcandidate_mass_4e = []
        # ZLLSSCRcandidate_pt_4e = []
        # ZLLSSCRcandidate_eta_4e = []
        # ZLLSSCRcandidate_phi_4e = []
        # ZLLSSCRcandidate_lep3_pt_4e = []
        # ZLLSSCRcandidate_lep3_eta_4e = []
        # ZLLSSCRcandidate_lep3_pdgId_4e = []
        # ZLLSSCRcandidate_lep4_pt_4e = []
        # ZLLSSCRcandidate_lep4_eta_4e = []
        # ZLLSSCRcandidate_lep4_pdgId_4e = []
        # ZLLSSCRcandidate_lep3_lostHits_4e = []
        # ZLLSSCRcandidate_lep4_lostHits_4e = []
        # ZLLSSCRcandidate_lep3_isGood_4e = []
        # ZLLSSCRcandidate_lep4_isGood_4e = []
        # ZLLSSCRcandidate_mass_4mu = []
        # ZLLSSCRcandidate_pt_4mu = []
        # ZLLSSCRcandidate_eta_4mu = []
        # ZLLSSCRcandidate_phi_4mu = []
        # ZLLSSCRcandidate_lep3_pt_4mu = []
        # ZLLSSCRcandidate_lep3_eta_4mu = []
        # ZLLSSCRcandidate_lep3_pdgId_4mu = []
        # ZLLSSCRcandidate_lep4_pt_4mu = []
        # ZLLSSCRcandidate_lep4_eta_4mu = []
        # ZLLSSCRcandidate_lep4_pdgId_4mu = []
        # ZLLSSCRcandidate_lep3_isGood_4mu = []
        # ZLLSSCRcandidate_lep4_isGood_4mu = []
        # ZLLSSCRcandidate_mass_2e2mu = []
        # ZLLSSCRcandidate_pt_2e2mu = []
        # ZLLSSCRcandidate_eta_2e2mu = []
        # ZLLSSCRcandidate_phi_2e2mu = []
        # ZLLSSCRcandidate_lep3_pt_2e2mu = []
        # ZLLSSCRcandidate_lep3_eta_2e2mu = []
        # ZLLSSCRcandidate_lep3_pdgId_2e2mu = []
        # ZLLSSCRcandidate_lep4_pt_2e2mu = []
        # ZLLSSCRcandidate_lep4_eta_2e2mu = []
        # ZLLSSCRcandidate_lep4_pdgId_2e2mu = []
        # ZLLSSCRcandidate_lep3_isGood_2e2mu = []
        # ZLLSSCRcandidate_lep4_isGood_2e2mu = []

        # for ZLLSSCRcandidate in event.ZLLSSCRcandidates:
        #     ZLLSSCRcandidate_mass.append(ZLLSSCRcandidate.mass)
        #     ZLLSSCRcandidate_pt.append(ZLLSSCRcandidate.pt)
        #     ZLLSSCRcandidate_eta.append(ZLLSSCRcandidate.eta)
        #     ZLLSSCRcandidate_phi.append(ZLLSSCRcandidate.phi)

        #     # Extract PDG IDs
        #     lep_ids = {
        #         abs(ZLLSSCRcandidate.Z1.lep1.pdgId),
        #         abs(ZLLSSCRcandidate.Z1.lep2.pdgId),
        #         abs(ZLLSSCRcandidate.Z2.lep1.pdgId),
        #         abs(ZLLSSCRcandidate.Z2.lep2.pdgId)
        #     }

        #     # Classify by decay channel
        #     if lep_ids == {11}:  
        #         ZLLSSCRcandidate_mass_4e.append(ZLLSSCRcandidate.mass)
        #         ZLLSSCRcandidate_pt_4e.append(ZLLSSCRcandidate.pt)
        #         ZLLSSCRcandidate_eta_4e.append(ZLLSSCRcandidate.eta)
        #         ZLLSSCRcandidate_phi_4e.append(ZLLSSCRcandidate.phi)
        #         ZLLSSCRcandidate_lep3_pt_4e.append(ZLLSSCRcandidate.Z2.lep1.pt)
        #         ZLLSSCRcandidate_lep3_eta_4e.append(ZLLSSCRcandidate.Z2.lep1.eta)
        #         ZLLSSCRcandidate_lep3_pdgId_4e.append(ZLLSSCRcandidate.Z2.lep1.pdgId)
        #         ZLLSSCRcandidate_lep4_pt_4e.append(ZLLSSCRcandidate.Z2.lep2.pt)
        #         ZLLSSCRcandidate_lep4_eta_4e.append(ZLLSSCRcandidate.Z2.lep2.eta)
        #         ZLLSSCRcandidate_lep4_pdgId_4e.append(ZLLSSCRcandidate.Z2.lep2.pdgId)
        #         ZLLSSCRcandidate_lep3_lostHits_4e.append(ZLLSSCRcandidate.Z2.lep1.lostHits)
        #         ZLLSSCRcandidate_lep4_lostHits_4e.append(ZLLSSCRcandidate.Z2.lep2.lostHits)
        #         ZLLSSCRcandidate_lep3_isGood_4e.append(ZLLSSCRcandidate.Z2.lep1.isFullID)
        #         ZLLSSCRcandidate_lep4_isGood_4e.append(ZLLSSCRcandidate.Z2.lep2.isFullID)
        #     elif lep_ids == {13}:  
        #         ZLLSSCRcandidate_mass_4mu.append(ZLLSSCRcandidate.mass)
        #         ZLLSSCRcandidate_pt_4mu.append(ZLLSSCRcandidate.pt)
        #         ZLLSSCRcandidate_eta_4mu.append(ZLLSSCRcandidate.eta)
        #         ZLLSSCRcandidate_phi_4mu.append(ZLLSSCRcandidate.phi)
        #         ZLLSSCRcandidate_lep3_pt_4mu.append(ZLLSSCRcandidate.Z2.lep1.pt)
        #         ZLLSSCRcandidate_lep3_eta_4mu.append(ZLLSSCRcandidate.Z2.lep1.eta)
        #         ZLLSSCRcandidate_lep3_pdgId_4mu.append(ZLLSSCRcandidate.Z2.lep1.pdgId)
        #         ZLLSSCRcandidate_lep4_pt_4mu.append(ZLLSSCRcandidate.Z2.lep2.pt)
        #         ZLLSSCRcandidate_lep4_eta_4mu.append(ZLLSSCRcandidate.Z2.lep2.eta)
        #         ZLLSSCRcandidate_lep4_pdgId_4mu.append(ZLLSSCRcandidate.Z2.lep2.pdgId)
        #         ZLLSSCRcandidate_lep3_isGood_4mu.append(ZLLSSCRcandidate.Z2.lep1.isFullID)
        #         ZLLSSCRcandidate_lep4_isGood_4mu.append(ZLLSSCRcandidate.Z2.lep2.isFullID)
        #     elif lep_ids == {11, 13}:  
        #         ZLLSSCRcandidate_mass_2e2mu.append(ZLLSSCRcandidate.mass)
        #         ZLLSSCRcandidate_pt_2e2mu.append(ZLLSSCRcandidate.pt)
        #         ZLLSSCRcandidate_eta_2e2mu.append(ZLLSSCRcandidate.eta)
        #         ZLLSSCRcandidate_phi_2e2mu.append(ZLLSSCRcandidate.phi)
        #         ZLLSSCRcandidate_lep3_pt_2e2mu.append(ZLLSSCRcandidate.Z2.lep1.pt)
        #         ZLLSSCRcandidate_lep3_eta_2e2mu.append(ZLLSSCRcandidate.Z2.lep1.eta)
        #         ZLLSSCRcandidate_lep3_pdgId_2e2mu.append(ZLLSSCRcandidate.Z2.lep1.pdgId)
        #         ZLLSSCRcandidate_lep4_pt_2e2mu.append(ZLLSSCRcandidate.Z2.lep2.pt)
        #         ZLLSSCRcandidate_lep4_eta_2e2mu.append(ZLLSSCRcandidate.Z2.lep2.eta)
        #         ZLLSSCRcandidate_lep4_pdgId_2e2mu.append(ZLLSSCRcandidate.Z2.lep2.pdgId)
        #         ZLLSSCRcandidate_lep3_isGood_2e2mu.append(ZLLSSCRcandidate.Z2.lep1.isFullID)
        #         ZLLSSCRcandidate_lep4_isGood_2e2mu.append(ZLLSSCRcandidate.Z2.lep2.isFullID)

        # # Store in output dictionary
        # out_data[self.ZLLSSCR_prefix + "mass"] = ZLLSSCRcandidate_mass
        # out_data[self.ZLLSSCR_prefix + "pt"] = ZLLSSCRcandidate_pt
        # out_data[self.ZLLSSCR_prefix + "eta"] = ZLLSSCRcandidate_eta 
        # out_data[self.ZLLSSCR_prefix + "phi"] = ZLLSSCRcandidate_phi

        # out_data[self.ZLLSSCR4e_prefix + "mass"] = ZLLSSCRcandidate_mass_4e
        # out_data[self.ZLLSSCR4e_prefix + "pt"] = ZLLSSCRcandidate_pt_4e
        # out_data[self.ZLLSSCR4e_prefix + "eta"] = ZLLSSCRcandidate_eta_4e
        # out_data[self.ZLLSSCR4e_prefix + "phi"] = ZLLSSCRcandidate_phi_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep3_pt"] = ZLLSSCRcandidate_lep3_pt_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep3_eta"] = ZLLSSCRcandidate_lep3_eta_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep3_pdgId"] = ZLLSSCRcandidate_lep3_pdgId_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep4_pt"] = ZLLSSCRcandidate_lep4_pt_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep4_eta"] = ZLLSSCRcandidate_lep4_eta_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep4_pdgId"] = ZLLSSCRcandidate_lep4_pdgId_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep3_lostHits"] = ZLLSSCRcandidate_lep3_lostHits_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep4_lostHits"] = ZLLSSCRcandidate_lep4_lostHits_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep3_isGood"] = ZLLSSCRcandidate_lep3_isGood_4e
        # out_data[self.ZLLSSCR4e_prefix + "lep4_isGood"] = ZLLSSCRcandidate_lep4_isGood_4e

        # out_data[self.ZLLSSCR4mu_prefix + "mass"] = ZLLSSCRcandidate_mass_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "pt"] = ZLLSSCRcandidate_pt_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "eta"] = ZLLSSCRcandidate_eta_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "phi"] = ZLLSSCRcandidate_phi_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep3_pt"] = ZLLSSCRcandidate_lep3_pt_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep3_eta"] = ZLLSSCRcandidate_lep3_eta_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep3_pdgId"] = ZLLSSCRcandidate_lep3_pdgId_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep4_pt"] = ZLLSSCRcandidate_lep4_pt_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep4_eta"] = ZLLSSCRcandidate_lep4_eta_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep4_pdgId"] = ZLLSSCRcandidate_lep4_pdgId_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep3_isGood"] = ZLLSSCRcandidate_lep3_isGood_4mu
        # out_data[self.ZLLSSCR4mu_prefix + "lep4_isGood"] = ZLLSSCRcandidate_lep4_isGood_4mu

        # out_data[self.ZLLSSCR2e2mu_prefix + "mass"] =  ZLLSSCRcandidate_mass_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "pt"] = ZLLSSCRcandidate_pt_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "eta"] = ZLLSSCRcandidate_eta_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "phi"] = ZLLSSCRcandidate_phi_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep3_pt"] = ZLLSSCRcandidate_lep3_pt_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep3_eta"] = ZLLSSCRcandidate_lep3_eta_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep3_pdgId"] = ZLLSSCRcandidate_lep3_pdgId_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep4_pt"] = ZLLSSCRcandidate_lep4_pt_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep4_eta"] = ZLLSSCRcandidate_lep4_eta_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep4_pdgId"] = ZLLSSCRcandidate_lep4_pdgId_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep3_isGood"] = ZLLSSCRcandidate_lep3_isGood_2e2mu
        # out_data[self.ZLLSSCR2e2mu_prefix + "lep4_isGood"] = ZLLSSCRcandidate_lep4_isGood_2e2mu


        # Z + L
        ZLcandidate_all_mass = []
        ZLcandidate_all_trimass = []
        ZLcandidate_all_pt = []
        ZLcandidate_all_eta = []
        ZLcandidate_all_phi = []
        ZLcandidate_all_pt2 = []
        ZLcandidate_all_eta2 = []
        ZLcandidate_all_sip3d = []
        ZLcandidate_all_dxy = []
        ZLcandidate_all_dz = []
        ZLcandidate_all_pfcand=[]
        ZLcandidate_all_pdgId=[]
        ZLcandidate_all_phi2=[]

        ZLcandidate_alle_mass = []
        ZLcandidate_alle_trimass = []
        ZLcandidate_alle_pt = []
        ZLcandidate_alle_eta = []
        ZLcandidate_alle_phi = []
        ZLcandidate_alle_pt2 = []
        ZLcandidate_alle_eta2 = []
        ZLcandidate_alle_sip3d = []
        ZLcandidate_alle_dxy = []
        ZLcandidate_alle_dz = []
        ZLcandidate_alle_pfcand=[]
        ZLcandidate_alle_pdgId=[]
        ZLcandidate_alle_phi2=[]
        ZLcandidate_alle_lostHits=[]


        ZLcandidate_allmu_mass = []
        ZLcandidate_allmu_trimass = []
        ZLcandidate_allmu_pt = []
        ZLcandidate_allmu_eta = []
        ZLcandidate_allmu_phi = []
        ZLcandidate_allmu_pt2 = []
        ZLcandidate_allmu_eta2 = []
        ZLcandidate_allmu_sip3d = []
        ZLcandidate_allmu_dxy = []
        ZLcandidate_allmu_dz = []
        ZLcandidate_allmu_pfcand=[]
        ZLcandidate_allmu_pdgId=[]
        ZLcandidate_allmu_phi2=[]


        for ZLcandidate_all in event.ZLcandidates_all:
            ZLcandidate_all_mass.append(ZLcandidate_all.mass)
            ZLcandidate_all_trimass.append(ZLcandidate_all.trimass)
            ZLcandidate_all_pt.append(ZLcandidate_all.pt)
            ZLcandidate_all_eta.append(ZLcandidate_all.eta)
            ZLcandidate_all_phi.append(ZLcandidate_all.phi)
            ZLcandidate_all_pt2.append(ZLcandidate_all.pt2)
            ZLcandidate_all_eta2.append(ZLcandidate_all.eta2)
            ZLcandidate_all_sip3d.append(ZLcandidate_all.sip3d)
            ZLcandidate_all_dz.append(ZLcandidate_all.dz)
            ZLcandidate_all_dxy.append(ZLcandidate_all.dxy)
            ZLcandidate_all_pfcand.append(ZLcandidate_all.pfcand)
            ZLcandidate_all_pdgId.append(ZLcandidate_all.pdgId)
            ZLcandidate_all_phi2.append(ZLcandidate_all.phi2)
        
        for ZLcandidate_alle in event.ZLcandidates_alle:
            ZLcandidate_alle_mass.append(ZLcandidate_alle.mass)
            ZLcandidate_alle_trimass.append(ZLcandidate_alle.trimass)
            ZLcandidate_alle_pt.append(ZLcandidate_alle.pt)
            ZLcandidate_alle_eta.append(ZLcandidate_alle.eta)
            ZLcandidate_alle_phi.append(ZLcandidate_alle.phi)
            ZLcandidate_alle_pt2.append(ZLcandidate_alle.pt2)
            ZLcandidate_alle_eta2.append(ZLcandidate_alle.eta2)
            ZLcandidate_alle_sip3d.append(ZLcandidate_alle.sip3d)
            ZLcandidate_alle_dz.append(ZLcandidate_alle.dz)
            ZLcandidate_alle_dxy.append(ZLcandidate_alle.dxy)
            ZLcandidate_alle_pfcand.append(ZLcandidate_alle.pfcand)  
            ZLcandidate_alle_pdgId.append(ZLcandidate_alle.pdgId)
            ZLcandidate_alle_phi2.append(ZLcandidate_alle.phi2)
            ZLcandidate_alle_lostHits.append(ZLcandidate_alle.lostHits)
        
        for ZLcandidate_allmu in event.ZLcandidates_allmu:
            ZLcandidate_allmu_mass.append(ZLcandidate_allmu.mass)
            ZLcandidate_allmu_trimass.append(ZLcandidate_allmu.trimass)
            ZLcandidate_allmu_pt.append(ZLcandidate_allmu.pt)
            ZLcandidate_allmu_eta.append(ZLcandidate_allmu.eta)
            ZLcandidate_allmu_phi.append(ZLcandidate_allmu.phi)
            ZLcandidate_allmu_pt2.append(ZLcandidate_allmu.pt2)
            ZLcandidate_allmu_eta2.append(ZLcandidate_allmu.eta2)
            ZLcandidate_allmu_sip3d.append(ZLcandidate_allmu.sip3d)
            ZLcandidate_allmu_dz.append(ZLcandidate_allmu.dz)
            ZLcandidate_allmu_dxy.append(ZLcandidate_allmu.dxy)
            ZLcandidate_allmu_pfcand.append(ZLcandidate_allmu.pfcand)
            ZLcandidate_allmu_pdgId.append(ZLcandidate_allmu.pdgId)
            ZLcandidate_allmu_phi2.append(ZLcandidate_allmu.phi2)
        

        out_data[self.ZLall_prefix + "mass"] = ZLcandidate_all_mass
        out_data[self.ZLall_prefix + "trimass"] = ZLcandidate_all_trimass
        out_data[self.ZLall_prefix + "pt"] = ZLcandidate_all_pt
        out_data[self.ZLall_prefix + "eta"] = ZLcandidate_all_eta 
        out_data[self.ZLall_prefix + "phi"] = ZLcandidate_all_phi
        out_data[self.ZLall_prefix + "pt2"] = ZLcandidate_all_pt2
        out_data[self.ZLall_prefix + "eta2"] = ZLcandidate_all_eta2
        out_data[self.ZLall_prefix + "sip3d"] = ZLcandidate_all_sip3d 
        out_data[self.ZLall_prefix + "dz"] = ZLcandidate_all_dz
        out_data[self.ZLall_prefix + "dxy"] = ZLcandidate_all_dxy
        out_data[self.ZLall_prefix + "pfcand"] = ZLcandidate_all_pfcand    
        out_data[self.ZLall_prefix + "pdgId"] = ZLcandidate_all_pdgId
        out_data[self.ZLall_prefix + "phi2"] = ZLcandidate_all_phi2



        out_data[self.ZLalle_prefix + "mass"] = ZLcandidate_alle_mass
        out_data[self.ZLalle_prefix + "trimass"] = ZLcandidate_alle_trimass
        out_data[self.ZLalle_prefix + "pt"] = ZLcandidate_alle_pt
        out_data[self.ZLalle_prefix + "eta"] = ZLcandidate_alle_eta 
        out_data[self.ZLalle_prefix + "phi"] = ZLcandidate_alle_phi
        out_data[self.ZLalle_prefix + "pt2"] = ZLcandidate_alle_pt2
        out_data[self.ZLalle_prefix + "eta2"] = ZLcandidate_alle_eta2
        out_data[self.ZLalle_prefix + "sip3d"] = ZLcandidate_alle_sip3d 
        out_data[self.ZLalle_prefix + "dz"] = ZLcandidate_alle_dz
        out_data[self.ZLalle_prefix + "dxy"] = ZLcandidate_alle_dxy
        out_data[self.ZLalle_prefix + "pfcand"] = ZLcandidate_alle_pfcand
        out_data[self.ZLalle_prefix + "pdgId"] = ZLcandidate_alle_pdgId
        out_data[self.ZLalle_prefix + "phi2"] = ZLcandidate_alle_phi2
        out_data[self.ZLalle_prefix + "lostHits"] = ZLcandidate_alle_lostHits


        out_data[self.ZLallmu_prefix + "mass"] = ZLcandidate_allmu_mass
        out_data[self.ZLallmu_prefix + "trimass"] = ZLcandidate_allmu_trimass
        out_data[self.ZLallmu_prefix + "pt"] = ZLcandidate_allmu_pt
        out_data[self.ZLallmu_prefix + "eta"] = ZLcandidate_allmu_eta 
        out_data[self.ZLallmu_prefix + "phi"] = ZLcandidate_allmu_phi
        out_data[self.ZLallmu_prefix + "pt2"] = ZLcandidate_allmu_pt2
        out_data[self.ZLallmu_prefix + "eta2"] = ZLcandidate_allmu_eta2
        out_data[self.ZLallmu_prefix + "sip3d"] = ZLcandidate_allmu_sip3d 
        out_data[self.ZLallmu_prefix + "dz"] = ZLcandidate_allmu_dz
        out_data[self.ZLallmu_prefix + "dxy"] = ZLcandidate_allmu_dxy
        out_data[self.ZLallmu_prefix + "pfcand"] = ZLcandidate_allmu_pfcand
        out_data[self.ZLallmu_prefix + "pdgId"] = ZLcandidate_allmu_pdgId
        out_data[self.ZLallmu_prefix + "phi2"] = ZLcandidate_allmu_phi2

        
        ZLcandidate_pass_mass = []
        ZLcandidate_pass_trimass = []
        ZLcandidate_pass_pt = []
        ZLcandidate_pass_eta = []
        ZLcandidate_pass_phi = []
        ZLcandidate_pass_pt2 = []
        ZLcandidate_pass_eta2 = []
        ZLcandidate_pass_sip3d = []
        ZLcandidate_pass_dxy = []
        ZLcandidate_pass_dz = []
        ZLcandidate_pass_pfcand=[]
        ZLcandidate_pass_pdgId=[]
        ZLcandidate_pass_phi2=[]

        ZLcandidate_passe_mass = []
        ZLcandidate_passe_trimass = []
        ZLcandidate_passe_pt = []
        ZLcandidate_passe_eta = []
        ZLcandidate_passe_phi = []
        ZLcandidate_passe_pt2 = []
        ZLcandidate_passe_eta2 = []
        ZLcandidate_passe_sip3d = []
        ZLcandidate_passe_dxy = []
        ZLcandidate_passe_dz = []
        ZLcandidate_passe_pfcand=[]
        ZLcandidate_passe_pdgId=[]
        ZLcandidate_passe_phi2=[]
        ZLcandidate_passe_lostHits=[]

        ZLcandidate_passmu_mass = []
        ZLcandidate_passmu_trimass = []
        ZLcandidate_passmu_pt = []
        ZLcandidate_passmu_eta = []
        ZLcandidate_passmu_phi = []
        ZLcandidate_passmu_pt2 = []
        ZLcandidate_passmu_eta2 = []
        ZLcandidate_passmu_sip3d = []
        ZLcandidate_passmu_dxy = []
        ZLcandidate_passmu_dz = []
        ZLcandidate_passmu_pfcand=[]
        ZLcandidate_passmu_pdgId=[]
        ZLcandidate_passmu_phi2=[]

        for ZLcandidate_pass in event.ZLcandidates_pass:
            ZLcandidate_pass_mass.append(ZLcandidate_pass.mass)
            ZLcandidate_pass_trimass.append(ZLcandidate_pass.trimass)
            ZLcandidate_pass_pt.append(ZLcandidate_pass.pt)
            ZLcandidate_pass_eta.append(ZLcandidate_pass.eta)
            ZLcandidate_pass_phi.append(ZLcandidate_pass.phi)
            ZLcandidate_pass_pt2.append(ZLcandidate_pass.pt2)
            ZLcandidate_pass_eta2.append(ZLcandidate_pass.eta2)
            ZLcandidate_pass_sip3d.append(ZLcandidate_pass.sip3d)
            ZLcandidate_pass_dz.append(ZLcandidate_pass.dz)
            ZLcandidate_pass_dxy.append(ZLcandidate_pass.dxy)
            ZLcandidate_pass_pfcand.append(ZLcandidate_pass.pfcand)
            ZLcandidate_pass_pdgId.append(ZLcandidate_pass.pdgId)
            ZLcandidate_pass_phi2.append(ZLcandidate_pass.phi2)
        
        for ZLcandidate_passe in event.ZLcandidates_passe:
            ZLcandidate_passe_mass.append(ZLcandidate_passe.mass)
            ZLcandidate_passe_trimass.append(ZLcandidate_passe.trimass)
            ZLcandidate_passe_pt.append(ZLcandidate_passe.pt)
            ZLcandidate_passe_eta.append(ZLcandidate_passe.eta)
            ZLcandidate_passe_phi.append(ZLcandidate_passe.phi)
            ZLcandidate_passe_pt2.append(ZLcandidate_passe.pt2)
            ZLcandidate_passe_eta2.append(ZLcandidate_passe.eta2)
            ZLcandidate_passe_sip3d.append(ZLcandidate_passe.sip3d)
            ZLcandidate_passe_dz.append(ZLcandidate_passe.dz)
            ZLcandidate_passe_dxy.append(ZLcandidate_passe.dxy)
            ZLcandidate_passe_pfcand.append(ZLcandidate_passe.pfcand)
            ZLcandidate_passe_pdgId.append(ZLcandidate_passe.pdgId)
            ZLcandidate_passe_phi2.append(ZLcandidate_passe.phi2)
            ZLcandidate_passe_lostHits.append(ZLcandidate_passe.lostHits)
        
        for ZLcandidate_passmu in event.ZLcandidates_passmu:
            ZLcandidate_passmu_mass.append(ZLcandidate_passmu.mass)
            ZLcandidate_passmu_trimass.append(ZLcandidate_passmu.trimass)
            ZLcandidate_passmu_pt.append(ZLcandidate_passmu.pt)
            ZLcandidate_passmu_eta.append(ZLcandidate_passmu.eta)
            ZLcandidate_passmu_phi.append(ZLcandidate_passmu.phi)
            ZLcandidate_passmu_pt2.append(ZLcandidate_passmu.pt2)
            ZLcandidate_passmu_eta2.append(ZLcandidate_passmu.eta2)
            ZLcandidate_passmu_sip3d.append(ZLcandidate_passmu.sip3d)
            ZLcandidate_passmu_dz.append(ZLcandidate_passmu.dz)
            ZLcandidate_passmu_dxy.append(ZLcandidate_passmu.dxy)
            ZLcandidate_passmu_pfcand.append(ZLcandidate_passmu.pfcand)
            ZLcandidate_passmu_pdgId.append(ZLcandidate_passmu.pdgId)
            ZLcandidate_passmu_phi2.append(ZLcandidate_passmu.phi2)
        

        out_data[self.ZLpass_prefix + "mass"] = ZLcandidate_pass_mass
        out_data[self.ZLpass_prefix + "trimass"] = ZLcandidate_pass_trimass
        out_data[self.ZLpass_prefix + "pt"] = ZLcandidate_pass_pt
        out_data[self.ZLpass_prefix + "eta"] = ZLcandidate_pass_eta 
        out_data[self.ZLpass_prefix + "phi"] = ZLcandidate_pass_phi
        out_data[self.ZLpass_prefix + "pt2"] = ZLcandidate_pass_pt2
        out_data[self.ZLpass_prefix + "eta2"] = ZLcandidate_pass_eta2
        out_data[self.ZLpass_prefix + "sip3d"] = ZLcandidate_pass_sip3d 
        out_data[self.ZLpass_prefix + "dz"] = ZLcandidate_pass_dz
        out_data[self.ZLpass_prefix + "dxy"] = ZLcandidate_pass_dxy
        out_data[self.ZLpass_prefix + "pfcand"] = ZLcandidate_pass_pfcand
        out_data[self.ZLpass_prefix + "pdgId"] = ZLcandidate_pass_pdgId
        out_data[self.ZLpass_prefix + "phi2"] = ZLcandidate_pass_phi2

        out_data[self.ZLpasse_prefix + "mass"] = ZLcandidate_passe_mass
        out_data[self.ZLpasse_prefix + "trimass"] = ZLcandidate_passe_trimass
        out_data[self.ZLpasse_prefix + "pt"] = ZLcandidate_passe_pt
        out_data[self.ZLpasse_prefix + "eta"] = ZLcandidate_passe_eta 
        out_data[self.ZLpasse_prefix + "phi"] = ZLcandidate_passe_phi
        out_data[self.ZLpasse_prefix + "pt2"] = ZLcandidate_passe_pt2
        out_data[self.ZLpasse_prefix + "eta2"] = ZLcandidate_passe_eta2
        out_data[self.ZLpasse_prefix + "sip3d"] = ZLcandidate_passe_sip3d 
        out_data[self.ZLpasse_prefix + "dz"] = ZLcandidate_passe_dz
        out_data[self.ZLpasse_prefix + "dxy"] = ZLcandidate_passe_dxy
        out_data[self.ZLpasse_prefix + "pfcand"] = ZLcandidate_passe_pfcand
        out_data[self.ZLpasse_prefix + "pdgId"] = ZLcandidate_passe_pdgId
        out_data[self.ZLpasse_prefix + "phi2"] = ZLcandidate_passe_phi2
        out_data[self.ZLpasse_prefix + "lostHits"] = ZLcandidate_passe_lostHits

        out_data[self.ZLpassmu_prefix + "mass"] = ZLcandidate_passmu_mass
        out_data[self.ZLpassmu_prefix + "trimass"] = ZLcandidate_passmu_trimass
        out_data[self.ZLpassmu_prefix + "pt"] = ZLcandidate_passmu_pt
        out_data[self.ZLpassmu_prefix + "eta"] = ZLcandidate_passmu_eta 
        out_data[self.ZLpassmu_prefix + "phi"] = ZLcandidate_passmu_phi
        out_data[self.ZLpassmu_prefix + "pt2"] = ZLcandidate_passmu_pt2
        out_data[self.ZLpassmu_prefix + "eta2"] = ZLcandidate_passmu_eta2
        out_data[self.ZLpassmu_prefix + "sip3d"] = ZLcandidate_passmu_sip3d 
        out_data[self.ZLpassmu_prefix + "dz"] = ZLcandidate_passmu_dz
        out_data[self.ZLpassmu_prefix + "dxy"] = ZLcandidate_passmu_dxy
        out_data[self.ZLpassmu_prefix + "pfcand"] = ZLcandidate_passmu_pfcand
        out_data[self.ZLpassmu_prefix + "pdgId"] = ZLcandidate_passmu_pdgId
        out_data[self.ZLpassmu_prefix + "phi2"] = ZLcandidate_passmu_phi2



        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True