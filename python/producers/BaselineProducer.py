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
    
    def __init__(self,lep1,lep2):
        self.lep1 = lep1
        self.lep2 = lep2
        self.pt = sumP4(self.lep1, self.lep2).Pt()
        self.eta = sumP4(self.lep1, self.lep2).Eta()
        self.phi = sumP4(self.lep1, self.lep2).Phi()
        self.mass = sumP4(self.lep1, self.lep2).M()

class BaselineProducer(Module):
    
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

        # Define the prefixes
        self.mu_prefix = "mu_"
        self.el_prefix = "el_"
        self.lep_prefix = "lep_"
        self.jet_prefix = "jet_"
        self.Z_prefix = "Z_"
        self.Zmu_prefix = "Zmu_"
        self.Zel_prefix = "Zel_"

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
            self.out.branch(self.Zmu_prefix + Z_var, "F", 20, lenVar="nZmu")
            self.out.branch(self.Zel_prefix + Z_var, "F", 20, lenVar="nZel")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if event.PV_npvsGood < 1: return False

        # Apply trigger selections 
        if self._select_triggers(event) is False:
            return False

        self._select_muons(event)
        self._select_electrons(event)
        self._select_jets(event)

        event.selectedLeptons = event.selectedMuons + event.selectedElectrons

        if not (
            (len(event.selectedMuons) == 2 and len(event.selectedElectrons) == 0) or
            (len(event.selectedElectrons) == 2 and len(event.selectedMuons) == 0)
        ):
            return False
        
        self._select_Z_candidates(event)

        if len(event.Zcandidates) > 1 or len(event.Zcandidates) == 0:
            return False

        self._fill_event_info(event)

        return True
        
    def _select_triggers(self, event):

        passTrigger = False 
        out_data = {}
        if self.year == "2022" or self.year == "2022EE" : # Checked that these are unprescaled in run 359751
            passSingleEle = event.HLT_Ele30_WPTight_Gsf #Note: we used Ele32 in 2018! 
            passSingleMu = event.HLT_IsoMu24
            # passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
            # passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            # passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            # passTriEle = False
            # passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        elif self.year == "2023" or self.year == "2023BPix" : # Checked that these are unprescaled, reference twikis for 2023 Eg & Muon Triggers https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary & https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2023
            passSingleEle = event.HLT_Ele30_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
            # passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
            # passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            # passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
            # passTriEle = False
            # passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        else:
            print(f"Year {self.year} not found")

        if self.isMC or self.sample == "any" :
            passTrigger = passSingleEle or passSingleMu # passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu
        else: # Data: ensure each event is taken only from a single sample
            if self.sample == "" : sys.exit("ERROR: sample must be set in data") # we may want to merge triggers for test runs 
            # if (self.sample in ["DoubleEle", "DoubleEG", "EGamma"] and (passDiEle or passTriEle)) or \
            # (self.sample in ["Muon", "DoubleMu", "DoubleMuon"] and (passDiMu or passTriMu) and not passDiEle and not passTriEle) or \
            # (self.sample in ["MuEG", "MuonEG"] and passMuEle and not passDiEle and not passTriEle and not passDiMu and not passTriMu) or \
            # (self.sample in ["SingleElectron", "EGamma"] and passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
            # (self.sample in ["SingleMuon", "Muon"] and passSingleMu and not passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle):
            #     passTrigger = True
            if (self.sample == "Muon" and passSingleMu) or \
               (self.sample == "EGamma" and passSingleEle):
                passTrigger = True

        # self.out.fillBranch("HLT_passZZ4lEle", passSingleEle or passDiEle or passTriEle)
        # self.out.fillBranch("HLT_passZZ4lMu", passSingleMu or passDiMu or passTriMu)
        # self.out.fillBranch("HLT_passZZ4lMuEle", passMuEle)
        self.out.fillBranch("HLT_passZZ4l", passTrigger)

        return passTrigger

    def _select_Z_candidates(self, event):

        event.Zcandidates = []

        """Select Z candidates from the event's selected muons."""
        event.Zcandidates_mu = []

        muon_pairs = itertools.combinations(event.selectedMuons, 2)

        for mu1, mu2 in muon_pairs:
            if mu1.pdgId == mu2.pdgId: # Select only opposite-sign muons
                continue

            dr_ll = deltaR(mu1.eta, mu1.phi, mu2.eta, mu2.phi)
            if dr_ll < 0.3:
                continue
            
            if not ((mu1.pt > 27 or mu2.pt > 27) and (mu1.pt > 15 or mu2.pt > 15)):
                continue

            Zcand_mu = Zcandidate(mu1, mu2)

            if Zcand_mu.mass < 60 or Zcand.mass > 120:
                continue
            
            event.Zcandidates_mu.append(Zcand_mu) 
            event.Zcandidates.append(Zcand_mu)

        """Select Z candidates from the event's selected electrons."""
        event.Zcandidates_el = []

        electron_pairs = itertools.combinations(event.selectedElectrons, 2)

        for el1, el2 in electron_pairs:
            if el1.pdgId == el2.pdgId: # Select only opposite-sign electrons
                continue

            dr_ll = deltaR(el1.eta, el1.phi, el2.eta, el2.phi)
            if dr_ll < 0.3:
                continue:
            
            if not ((el1.pt > 33 or el2.pt > 33) and (el1.pt > 15 or el2.pt > 15)):
                continue

            Zcand_el = Zcandidate(el1, el2)

            if Zcand_el.mass < 60 or Zcand.mass > 120:
                continue

            event.Zcandidates_el.append(Zcand_el)
            event.Zcandidates.append(Zcand_el)

    def _select_muons(self, event):

        event.selectedMuons = []

        muons = Collection(event, "Muon")
        
        for mu in muons:
            
            passMuID = mu.isPFcand or (mu.highPtId>0 and mu.pt>200)
            
            if mu.pt > 10 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1 mu.pfRelIso03_all < 0.15:
                mu._wp_ID = 'TightID'
                mu._wp_Iso = 'LooseRelIso'
                event.selectedMuons.append(mu)
            
     def _select_electrons(self, event):

        event.selectedElectrons = []

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 10 and abs(el.eta) < 2.4 and abs(el.dxy) < 0.5 and abs(el.dz) < 1:
                # el._wp_ID = 'wp90iso'
                
                # mva = el.mvaHZZIso
                # if etaSC < 0.8:
                #     if el.pt < 10 and mva < 0.9044286167: continue
                #     if el.pt >= 10 and mva < 0.1968600840: continue
                # elif 0.8 < etaSC < 1.479:
                #     if el.pt < 10 and mva < 0.9094166886: continue
                #     if el.pt >= 10 and mva < 0.0759172100: continue
                # else:
                #     if el.pt < 10 and mva < 0.9443653660: continue
                #     if el.pt >= 10 and mva < -0.5169136775: continue      
                                    
                event.selectedElectrons.append(el)

    def _select_jets(self, event):

        event.selectedJets = []

        jets = Collection(event, "Jet")
        FsrPhotons = Collection(event, "FsrPhoton")

        for jet in jets:
            if jet.pt <= 20 and abs(jet.eta) >= 2.5:
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
        ak4_pt = []
        ak4_eta = []
        ak4_phi = []
        ak4_mass = []
        ak4_hadronFlavour = []
        
        for jet in event.selectedJets:
            ak4_pt.append(jet.pt)
            ak4_eta.append(jet.eta)
            ak4_phi.append(jet.phi)
            ak4_mass.append(jet.mass)
            if self.isMC: ak4_hadronFlavour.append(jet.hadronFlavour)
        
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
                
        out_data[self.Z_prefix + "mass"] = Zcandidate_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_phi 

        ## Z(mu,mu) candidates
        Zcandidate_mu_mass = []
        Zcandidate_mu_pt = []
        Zcandidate_mu_eta = []
        Zcandidate_mu_phi = []
        for Zcandidate in event.Zcandidates_mu:
            Zcandidate_mu_mass.append(Zcandidate.mass)
            Zcandidate_mu_pt.append(Zcandidate.pt)
            Zcandidate_mu_eta.append(Zcandidate.eta)
            Zcandidate_mu_phi.append(Zcandidate.phi)
                
        out_data[self.Z_prefix + "mass"] = Zcandidate_mu_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_mu_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_mu_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_mu_phi

        ## Z(el,el) candidates
        Zcandidate_el_mass = []
        Zcandidate_el_pt = []
        Zcandidate_el_eta = []
        Zcandidate_el_phi = []
        for Zcandidate in event.Zcandidates_el:
            Zcandidate_el_mass.append(Zcandidate.mass)
            Zcandidate_el_pt.append(Zcandidate.pt)
            Zcandidate_el_eta.append(Zcandidate.eta)
            Zcandidate_el_phi.append(Zcandidate.phi)
                
        out_data[self.Z_prefix + "mass"] = Zcandidate_el_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_el_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_el_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_el_phi 

        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True