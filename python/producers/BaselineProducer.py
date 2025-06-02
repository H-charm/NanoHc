from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import ROOT
import math
import itertools
from functools import cmp_to_key
from ..helpers.utils import sumP4
ROOT.PyConfig.IgnoreCommandLineOptions = True

lumi_dict = {"2022": 7.98, "2022EE": 26.67, "2023": 17.794, "2023BPix": 9.451}

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
        self.mass2 = sumP4(self.Z1, self.Z2).M()

class BaselineProducer(Module):
    
    def __init__(self, year, dataset_type, sample):
        self.year = year
        self.sample = sample
        self.dataset_type = dataset_type

        # Define the variables you want to plot
        self.lep_vars = ["pt","eta","phi", "pdgId"]
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

        # Define the prefixes
        self.mu_prefix = "mu_"
        self.el_prefix = "el_"
        self.lep_prefix = "lep_"
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

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        
        self.out = wrappedOutputTree
        
        # Define lepton branches
        for lep_var in self.lep_vars:
            self.out.branch(self.mu_prefix + lep_var, "F", 20, lenVar="nMu")
            self.out.branch(self.el_prefix + lep_var, "F", 20, lenVar="nEl")
            self.out.branch(self.lep_prefix + lep_var, "F", 20, lenVar="nLep")
        
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

        # # apply trigger selections on data
        # if not self.isMC: 
        #     if self._select_triggers(event) is False:
        #         return False

        # Apply trigger selections on trigger and data
        if self._select_triggers(event) is False:
            return False
        
        self._associate_fsr_photons_and_leptons(event)
        # self._select_muons(event)
        # self._select_electrons(event)  
        event.selectedLeptons = event.selectedMuons + event.selectedElectrons
        if len(event.selectedLeptons) < 4: return False

        self._select_jets(event)
        # if len(event.selectedJets) == 0:
        #     return False
        
        self._select_Z_candidates(event)
        if len(event.Zcandidates) < 2:
            return False
        
        self._select_ZZ_candidates(event)
        
        self._select_H_candidates(event)
        
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

        # if not out_data['passTriggers']:
        #     return False

        # for key in out_data:
        #     self.out.fillBranch(key, out_data[key]

        self.out.fillBranch("HLT_passZZ4lEle", passSingleEle or passDiEle or passTriEle)
        self.out.fillBranch("HLT_passZZ4lMu", passSingleMu or passDiMu or passTriMu)
        self.out.fillBranch("HLT_passZZ4lMuEle", passMuEle)
        self.out.fillBranch("HLT_passZZ4l", passTrigger)

        return passTrigger

    def _select_Z_candidates(self, event):

        event.Zcandidates = []
        
        # mva_leptons = [lepton for lepton in event.selectedLeptons if lepton.mvaTOP > 0.9]
        # lepton_pairs = list(itertools.combinations(mva_leptons, 2))
        lepton_pairs = list(itertools.combinations(event.selectedLeptons, 2))

        
        for lepton_pair in lepton_pairs:
            
            # we need same flavor and opposite charge
            if (lepton_pair[0].pdgId + lepton_pair[1].pdgId) != 0:
                continue
            
            # let's put negative charged lepton always as lep1 (useful later)
            lep1, lep2 = (lepton_pair[0], lepton_pair[1]) if lepton_pair[0].pdgId < 0 else (lepton_pair[1], lepton_pair[0])

            fsrPhotons = Collection(event, "FsrPhoton")
            fsrIndices = {
                "muFsrPhotonIdx": getattr(event, "muFsrPhotonIdx", []),
                "eleFsrPhotonIdx": getattr(event, "eleFsrPhotonIdx", [])
            }
            Zcand = Zcandidate(lep1, lep2, fsrPhotons, fsrIndices)
            
            if Zcand.mass < 12 or Zcand.mass > 120:
                continue
        
            event.Zcandidates.append(Zcand)

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
                
        Zcand_pairs = list(itertools.combinations(event.Zcandidates, 2))
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

            # Trigger pT cut
            leptons = [Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2]
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
                    return 1   # b is better

        event.Hcandidates = []

        if len(event.ZZcandidates) == 0:
            return
        ##---------------------
        ## Comparator
        best_candidate = min(event.ZZcandidates, key=cmp_to_key(best_candidate_comparator))
        event.Hcandidates.append(best_candidate)


        # mZ = 91.1876
    
        # # Filter ZZ candidates within Higgs mass window
        # valid_ZZ_candidates = []
        # for ZZcand in event.ZZcandidates:
        #     # if 105 <= ZZcand.mass <= 140:
        #     valid_ZZ_candidates.append(ZZcand)
        
        # if not valid_ZZ_candidates:
        #     return  # No valid Higgs candidate found
        # else:
        #     if 
        #     best_candidate = min(
        #         valid_ZZ_candidates,
        #         key=lambda ZZcand: abs(ZZcand.Z1.mass - mZ) 
        #     )

        # event.Hcandidates.append(best_candidate)
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
    # def _select_muons(self, event):
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        # muFsrPhotonIdx = getattr(event, "muFsrPhotonIdx", [-1] * len(muons))

        # selectedFSR = [fsrPhotons[i] for i in muFsrPhotonIdx if i >= 0]

        event.selectedMuons = []

        for imu, mu in enumerate(muons):

            isoCorr = self.isoFsrCorr(mu, selectedFSR)
            passIso = isoCorr < 0.35

            # full ID
            passID = (mu.isPFcand or (mu.highPtId > 0 and mu.pt > 200))
            passKinematics = (
                mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5
                and abs(mu.dz) < 1 and abs(mu.sip3d) < 4
            )

            if passID and passKinematics and passIso and (mu.isGlobal or (mu.isTracker and mu.nStations > 0)):
                mu._wp_Iso = 'LoosePFIso'
                mu._wp_ID = 'TightID'
                mu.iso = isoCorr
                mu.index = imu
                event.selectedMuons.append(mu)

            
    #     ## taken from here https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/nanoZZ4lAnalysis.py       
    # def _select_electrons(self, event):
        electrons = Collection(event, "Electron")
        fsrPhotons = Collection(event, "FsrPhoton")
        # eleFsrPhotonIdx = getattr(event, "eleFsrPhotonIdx", [-1] * len(electrons))

        # selectedFSR = [fsrPhotons[i] for i in eleFsrPhotonIdx if i >= 0]

        event.selectedElectrons = []

        for iel, el in enumerate(electrons):
            etaSC = el.eta + el.deltaEtaSC
            isoCorr = self.isoFsrCorr(el, selectedFSR)

            if el.pt > 7 and abs(el.eta) < 2.5 and abs(el.dxy) < 0.5 and abs(el.dz) < 1 and abs(el.sip3d) < 4:
                # full BDT ID
                el.etaSC = etaSC
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
                
                el._wp_ID = 'wp90iso'
                el.iso = isoCorr
                el.etaSC = etaSC
                el.index = iel
                event.selectedElectrons.append(el)



            # baseline cut

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
                    
        ## jets  
        ak4_bdisc = []
        ak4_cvbdisc = []
        ak4_cvldisc = []
        ak4_gvudsdisc = []
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
        for Zcandidate in event.Zcandidates:
            Zcandidate_mass.append(Zcandidate.mass)
            Zcandidate_pt.append(Zcandidate.pt)
            Zcandidate_eta.append(Zcandidate.eta)
            Zcandidate_phi.append(Zcandidate.phi)
            if Zcandidate.is_onshell: Zcandidate_onshell_mass.append(Zcandidate.mass)
            else: Zcandidate_offshell_mass.append(Zcandidate.mass)

                
        out_data[self.Z_prefix + "mass"] = Zcandidate_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_phi 
        out_data[self.Z_prefix + "onshell_mass"] = Zcandidate_onshell_mass
        out_data[self.Z_prefix + "offshell_mass"] = Zcandidate_offshell_mass

        ## ZZ candidates
        # Initialize output lists
        ZZcandidate_mass = []
        ZZcandidate_mass_4e=[]
        ZZcandidate_mass_4mu=[]
        ZZcandidate_mass_2e2mu=[]
        ZZcandidate_pt = []
        ZZcandidate_eta = []
        ZZcandidate_phi = []

        for ZZcandidate in event.ZZcandidates:
            ZZcandidate_mass.append(ZZcandidate.mass)
            ZZcandidate_pt.append(ZZcandidate.pt)
            ZZcandidate_eta.append(ZZcandidate.eta)
            ZZcandidate_phi.append(ZZcandidate.phi)

            # Extract PDG IDs
            lep_ids = {
                abs(ZZcandidate.Z1.lep1.pdgId),
                abs(ZZcandidate.Z1.lep2.pdgId),
                abs(ZZcandidate.Z2.lep1.pdgId),
                abs(ZZcandidate.Z2.lep2.pdgId)
            }

            # Classify by decay channel
            if lep_ids == {11}:  
                ZZcandidate_mass_4e.append(ZZcandidate.mass)
            elif lep_ids == {13}:  
                ZZcandidate_mass_4mu.append(ZZcandidate.mass)
            elif lep_ids == {11, 13}:  
                ZZcandidate_mass_2e2mu.append(ZZcandidate.mass)

        # Store in output dictionary
        out_data[self.ZZ_prefix + "mass"] = ZZcandidate_mass
        out_data[self.ZZ4e_prefix + "mass"] = ZZcandidate_mass_4e
        out_data[self.ZZ4mu_prefix + "mass"] = ZZcandidate_mass_4mu
        out_data[self.ZZ2e2mu_prefix + "mass"] = ZZcandidate_mass_2e2mu
        out_data[self.ZZ_prefix + "pt"] = ZZcandidate_pt
        out_data[self.ZZ_prefix + "eta"] = ZZcandidate_eta
        out_data[self.ZZ_prefix + "phi"] = ZZcandidate_phi

        # Similar structure for Higgs candidates
        Hcandidate_mass = []
        Hcandidate_mass_4e=[]
        Hcandidate_mass_4mu=[]
        Hcandidate_mass_2e2mu=[]
        Hcandidate_pt = []
        Hcandidate_eta = []
        Hcandidate_phi = []

        for Hcandidate in event.Hcandidates:
            Hcandidate_mass.append(Hcandidate.mass)
            Hcandidate_pt.append(Hcandidate.pt)
            Hcandidate_eta.append(Hcandidate.eta)
            Hcandidate_phi.append(Hcandidate.phi)

            # Extract PDG IDs
            lep_ids = {
                abs(Hcandidate.Z1.lep1.pdgId),
                abs(Hcandidate.Z1.lep2.pdgId),
                abs(Hcandidate.Z2.lep1.pdgId),
                abs(Hcandidate.Z2.lep2.pdgId)
            }

            # Classify by decay channel
            if lep_ids == {11}:  
                Hcandidate_mass_4e.append(Hcandidate.mass)
            elif lep_ids == {13}:  
                Hcandidate_mass_4mu.append(Hcandidate.mass)
            elif lep_ids == {11, 13}:  
                Hcandidate_mass_2e2mu.append(Hcandidate.mass)
                
        out_data[self.H_prefix + "mass"] = Hcandidate_mass
        out_data[self.H4e_prefix + "mass"] = Hcandidate_mass_4e
        out_data[self.H4mu_prefix + "mass"] = Hcandidate_mass_4mu
        out_data[self.H2e2mu_prefix + "mass"] =  Hcandidate_mass_2e2mu
        out_data[self.H_prefix + "pt"] = Hcandidate_pt
        out_data[self.H_prefix + "eta"] = Hcandidate_eta 
        out_data[self.H_prefix + "phi"] = Hcandidate_phi 
                
        # if self.isMC:
        #     out_data["l1PreFiringWeight"] = event.L1PreFiringWeight_Nom                
                
        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True