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

        self.lep1 = lep1
        self.lep2 = lep2
        self.pt_2 = sumP4(self.lep1, self.lep2).Pt()
        self.eta_2 = sumP4(self.lep1, self.lep2).Eta()
        self.phi_2 = sumP4(self.lep1, self.lep2).Phi()
        self.mass_2 = sumP4(self.lep1, self.lep2).M()
        self.dR = deltaR(lep1.eta, lep1.phi, lep2.eta, lep2.phi)
        self.deta = abs(lep1.eta - lep2.eta)
        self.dphi = abs(lep1.phi - lep2.phi)

        self.lep1_pt = lep1.pt
        self.lep1_eta = lep1.eta
        self.lep1_phi = lep1.phi

        self.lep2_pt = lep2.pt
        self.lep2_eta = lep2.eta
        self.lep2_phi = lep2.phi

    def sumpt(self):
        return self.lep1.pt + self.lep2.pt

    def finalState(self):
        return self.lep1.pdgId * self.lep2.pdgId


class BaselineProducer(Module):
    
    def __init__(self, year, dataset_type, sample):
        self.year = year
        self.sample = sample
        self.dataset_type = dataset_type

        # Define the variables you want to plot
        self.lep_vars = ["pt","eta","phi","pdgId"]
        self.Z_vars = ["pt","eta","phi","mass","dR", "deta", "dphi","lep1_pt", "lep1_eta", "lep1_phi", "lep2_pt", "lep2_eta", "lep2_phi", "mass_2", "phi_2", "eta_2", "pt_2"]
        self.jet_vars = ["pt", "eta", "phi", "number"]

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
        
        for jet_var in self.jet_vars:
            self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar = "nJet")
        
        # Define trigger branches
        self.out.branch("HLT_pass", "O")      # pass trigger requirements for the given sample (including sample precedence vetos) 

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

        self._associate_fsr_photons_and_leptons(event)

        # self._select_muons(event)
        # self._select_electrons(event)
        event.selectedLeptons = event.selectedMuons + event.selectedElectrons
        self._select_jets(event)

        if not (
            (len(event.selectedMuons) == 2 and len(event.selectedElectrons) == 0) or
            (len(event.selectedElectrons) == 2 and len(event.selectedMuons) == 0)
        ):
            return False
        
        self._select_Z_candidates(event)

        if len(event.Zcandidates) != 1:
            return False

        self._fill_event_info(event)

        return True
        
    def _select_triggers(self, event):

        passTrigger = False 
        out_data = {}
        if self.year == "2022" or self.year == "2022EE" : # Checked that these are unprescaled in run 359751
            passSingleEle = event.HLT_Ele30_WPTight_Gsf 
            passSingleMu = event.HLT_IsoMu24
        elif self.year == "2023" or self.year == "2023BPix" : # Checked that these are unprescaled, reference twikis for 2023 Eg & Muon Triggers https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary & https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2023
            passSingleEle = event.HLT_Ele30_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
        else:
            print(f"Year {self.year} not found")

        if self.isMC or self.sample == "any" :
            passTrigger = passSingleEle or passSingleMu # passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu
        else: # Data: ensure each event is taken only from a single sample
            if self.sample == "" : sys.exit("ERROR: sample must be set in data") # we may want to merge triggers for test runs 
            if (self.sample == "Muon" and passSingleMu) or \
               (self.sample == "SingleMuon" and passSingleMu) or \
               (self.sample == "EGamma" and passSingleEle):
                passTrigger = True

        self.out.fillBranch("HLT_pass", passTrigger)

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
            
            if not ((mu1.pt > 27 and mu2.pt > 15) or (mu1.pt > 15 and mu2.pt > 27)):
                continue

            fsrPhotons = Collection(event, "FsrPhoton")
            fsrIndices = {
                "muFsrPhotonIdx": getattr(event, "muFsrPhotonIdx", [])
            }

            if mu1.pt > mu2.pt:
                Zcand_mu = Zcandidate(mu1, mu2, fsrPhotons, fsrIndices)
            else:
                Zcand_mu = Zcandidate(mu2, mu1, fsrPhotons, fsrIndices)

            if Zcand_mu.mass < 60 or Zcand_mu.mass > 120:
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
                continue
            
            if not ((el1.pt > 33 and el2.pt > 15) or (el1.pt > 15 and el2.pt > 33)):
                continue
            
            fsrPhotons = Collection(event, "FsrPhoton")
            fsrIndices = {
                "eleFsrPhotonIdx": getattr(event, "eleFsrPhotonIdx", [])
            }

            if el1.pt > el2.pt:
                Zcand_el = Zcandidate(el1, el2, fsrPhotons, fsrIndices)
            else:
                Zcand_el = Zcandidate(el2, el1, fsrPhotons, fsrIndices)

            if Zcand_el.mass < 60 or Zcand_el.mass > 120:
                continue

            event.Zcandidates_el.append(Zcand_el)
            event.Zcandidates.append(Zcand_el)

    # def _select_muons(self, event):

    #     event.selectedMuons = []

    #     muons = Collection(event, "Muon")
        
    #     for mu in muons:
            
    #         passMuID = mu.isPFcand or (mu.highPtId>0 and mu.pt>200)
            
    #         if mu.pt > 10 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1 and mu.pfRelIso03_all < 0.15 and mu.tightId == True:
    #             mu._wp_ID = 'TightID'
    #             mu._wp_Iso = 'TightPFIso'
    #             event.selectedMuons.append(mu)
            
    # def _select_electrons(self, event):
    #     event.selectedElectrons = []

    #     electrons = Collection(event, "Electron")
    #     for el in electrons:
    #         el.etaSC = el.eta + el.deltaEtaSC
    #         if el.pt > 10 and abs(el.eta) < 2.4 and abs(el.dxy) < 0.5 and abs(el.dz) < 1 and el.mvaIso_WP80 == True:
    #             el._wp_ID = 'wp80iso'
    #             event.selectedElectrons.append(el)
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
                
                el.iso = isoCorr
                el._wp_ID = 'wp90iso'
                el.etaSC = etaSC
                el.index = iel
                event.selectedElectrons.append(el)



    def _select_jets(self, event):

        event.selectedJets = []

        jets = Collection(event, "Jet")
        FsrPhotons = Collection(event, "FsrPhoton")

        for jet in jets:
            if jet.pt <= 20 and abs(jet.eta) >= 2.5 and jet.jetId != 6:
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

        ## Jets
        jet_pt = []
        jet_eta = []
        jet_phi = []
        jet_n = []

        for jet in event.selectedJets:
            jet_pt.append(jet.pt)
            jet_eta.append(jet.eta)
            jet_phi.append(jet.phi)
            #jet_n.append(nJet)
        out_data[self.jet_prefix + "pt"] = jet_pt
        out_data[self.jet_prefix + "eta"] = jet_eta
        out_data[self.jet_prefix + "phi"] = jet_phi
        #out_data[self.jet_prefix + "number"] = jet_n
                   
        ## Z candidates
        Zcandidate_mass = []
        Zcandidate_pt = []
        Zcandidate_eta = []
        Zcandidate_phi = []
        Zcandidate_dR = []
        Zcandidate_deta = []
        Zcandidate_dphi = []
        Zcandidate_lep1_pt = []
        Zcandidate_lep1_eta = []
        Zcandidate_lep1_phi = []
        Zcandidate_lep2_pt = []
        Zcandidate_lep2_eta = []
        Zcandidate_lep2_phi = []
        Zcandidate_mass_2 = []
        Zcandidate_pt_2 = []
        Zcandidate_eta_2 = []
        Zcandidate_phi_2 = []

        for Zcandidate in event.Zcandidates:
            Zcandidate_mass.append(Zcandidate.mass)
            Zcandidate_pt.append(Zcandidate.pt)
            Zcandidate_eta.append(Zcandidate.eta)
            Zcandidate_phi.append(Zcandidate.phi)
            Zcandidate_dR.append(Zcandidate.dR)
            Zcandidate_deta.append(Zcandidate.deta)
            Zcandidate_dphi.append(Zcandidate.dphi)
            Zcandidate_lep1_pt.append(Zcandidate.lep1_pt)
            Zcandidate_lep1_eta.append(Zcandidate.lep1_eta)
            Zcandidate_lep1_phi.append(Zcandidate.lep1_phi)
            Zcandidate_lep2_pt.append(Zcandidate.lep2_pt)
            Zcandidate_lep2_eta.append(Zcandidate.lep2_eta)
            Zcandidate_lep2_phi.append(Zcandidate.lep2_phi)
            Zcandidate_mass_2.append(Zcandidate.mass_2)
            Zcandidate_pt_2.append(Zcandidate.pt_2)
            Zcandidate_eta_2.append(Zcandidate.eta_2)
            Zcandidate_phi_2.append(Zcandidate.phi_2)

        out_data[self.Z_prefix + "mass"] = Zcandidate_mass
        out_data[self.Z_prefix + "pt"] = Zcandidate_pt
        out_data[self.Z_prefix + "eta"] = Zcandidate_eta 
        out_data[self.Z_prefix + "phi"] = Zcandidate_phi
        out_data[self.Z_prefix + "dR"] = Zcandidate_dR
        out_data[self.Z_prefix + "deta"] = Zcandidate_deta
        out_data[self.Z_prefix + "dphi"] = Zcandidate_dphi
        out_data[self.Z_prefix + "lep1_pt"] =  Zcandidate_lep1_pt
        out_data[self.Z_prefix + "lep1_eta"] =  Zcandidate_lep1_eta
        out_data[self.Z_prefix + "lep1_phi"] =  Zcandidate_lep1_phi
        out_data[self.Z_prefix + "lep2_pt"] =  Zcandidate_lep2_pt
        out_data[self.Z_prefix + "lep2_eta"] =  Zcandidate_lep2_eta
        out_data[self.Z_prefix + "lep2_phi"] =  Zcandidate_lep2_phi
        out_data[self.Z_prefix + "mass_2"] = Zcandidate_mass_2
        out_data[self.Z_prefix + "pt_2"] = Zcandidate_pt_2
        out_data[self.Z_prefix + "eta_2"] = Zcandidate_eta_2 
        out_data[self.Z_prefix + "phi_2"] = Zcandidate_phi_2

        ## Z(mu,mu) candidates
        Zcandidate_mu_mass = []
        Zcandidate_mu_pt = []
        Zcandidate_mu_eta = []
        Zcandidate_mu_phi = []
        Zcandidate_mu_dR = []
        Zcandidate_mu_deta = []
        Zcandidate_mu_dphi = []
        Zcandidate_mu_lep1_pt = []
        Zcandidate_mu_lep1_eta = []
        Zcandidate_mu_lep1_phi = []
        Zcandidate_mu_lep2_pt = []
        Zcandidate_mu_lep2_eta = []
        Zcandidate_mu_lep2_phi = []
        Zcandidate_mu_mass_2 = []
        Zcandidate_mu_pt_2 = []
        Zcandidate_mu_eta_2 = []
        Zcandidate_mu_phi_2 = []
        for Zcandidate in event.Zcandidates_mu:
            Zcandidate_mu_mass.append(Zcandidate.mass)
            Zcandidate_mu_pt.append(Zcandidate.pt)
            Zcandidate_mu_eta.append(Zcandidate.eta)
            Zcandidate_mu_phi.append(Zcandidate.phi)
            Zcandidate_mu_dR.append(Zcandidate.dR)
            Zcandidate_mu_deta.append(Zcandidate.deta)
            Zcandidate_mu_dphi.append(Zcandidate.dphi)
            Zcandidate_mu_lep1_pt.append(Zcandidate.lep1_pt)
            Zcandidate_mu_lep1_eta.append(Zcandidate.lep1_eta)
            Zcandidate_mu_lep1_phi.append(Zcandidate.lep1_phi)
            Zcandidate_mu_lep2_pt.append(Zcandidate.lep2_pt)
            Zcandidate_mu_lep2_eta.append(Zcandidate.lep2_eta)
            Zcandidate_mu_lep2_phi.append(Zcandidate.lep2_phi)
            Zcandidate_mu_mass_2.append(Zcandidate.mass_2)
            Zcandidate_mu_pt_2.append(Zcandidate.pt_2)
            Zcandidate_mu_eta_2.append(Zcandidate.eta_2)
            Zcandidate_mu_phi_2.append(Zcandidate.phi_2)
                
        out_data[self.Zmu_prefix + "mass"] = Zcandidate_mu_mass
        out_data[self.Zmu_prefix + "pt"] = Zcandidate_mu_pt
        out_data[self.Zmu_prefix + "eta"] = Zcandidate_mu_eta 
        out_data[self.Zmu_prefix + "phi"] = Zcandidate_mu_phi
        out_data[self.Zmu_prefix + "dR"] = Zcandidate_mu_dR
        out_data[self.Zmu_prefix + "deta"] = Zcandidate_mu_deta
        out_data[self.Zmu_prefix + "dphi"] = Zcandidate_mu_dphi
        out_data[self.Zmu_prefix + "lep1_pt"] =  Zcandidate_mu_lep1_pt
        out_data[self.Zmu_prefix + "lep1_eta"] =  Zcandidate_mu_lep1_eta
        out_data[self.Zmu_prefix + "lep1_phi"] =  Zcandidate_mu_lep1_phi
        out_data[self.Zmu_prefix + "lep2_pt"] =  Zcandidate_mu_lep2_pt
        out_data[self.Zmu_prefix + "lep2_eta"] =  Zcandidate_mu_lep2_eta
        out_data[self.Zmu_prefix + "lep2_phi"] =  Zcandidate_mu_lep2_phi
        out_data[self.Zmu_prefix + "mass_2"] = Zcandidate_mu_mass_2
        out_data[self.Zmu_prefix + "pt_2"] = Zcandidate_mu_pt_2
        out_data[self.Zmu_prefix + "eta_2"] = Zcandidate_mu_eta_2 
        out_data[self.Zmu_prefix + "phi_2"] = Zcandidate_mu_phi_2

        ## Z(el,el) candidates
        Zcandidate_el_mass = []
        Zcandidate_el_pt = []
        Zcandidate_el_eta = []
        Zcandidate_el_phi = []
        Zcandidate_el_dR = []
        Zcandidate_el_deta = []
        Zcandidate_el_dphi = []
        Zcandidate_el_lep1_pt = []
        Zcandidate_el_lep1_eta = []
        Zcandidate_el_lep1_phi = []
        Zcandidate_el_lep2_pt = []
        Zcandidate_el_lep2_eta = []
        Zcandidate_el_lep2_phi = []
        Zcandidate_el_mass_2 = []
        Zcandidate_el_pt_2 = []
        Zcandidate_el_eta_2 = []
        Zcandidate_el_phi_2 = []
        for Zcandidate in event.Zcandidates_el:
            Zcandidate_el_mass.append(Zcandidate.mass)
            Zcandidate_el_pt.append(Zcandidate.pt)
            Zcandidate_el_eta.append(Zcandidate.eta)
            Zcandidate_el_phi.append(Zcandidate.phi)
            Zcandidate_el_dR.append(Zcandidate.dR)
            Zcandidate_el_deta.append(Zcandidate.deta)
            Zcandidate_el_dphi.append(Zcandidate.dphi)
            Zcandidate_el_lep1_pt.append(Zcandidate.lep1_pt)
            Zcandidate_el_lep1_eta.append(Zcandidate.lep1_eta)
            Zcandidate_el_lep1_phi.append(Zcandidate.lep1_phi)
            Zcandidate_el_lep2_pt.append(Zcandidate.lep2_pt)
            Zcandidate_el_lep2_eta.append(Zcandidate.lep2_eta)
            Zcandidate_el_lep2_phi.append(Zcandidate.lep2_phi)
            Zcandidate_el_mass_2.append(Zcandidate.mass_2)
            Zcandidate_el_pt_2.append(Zcandidate.pt_2)
            Zcandidate_el_eta_2.append(Zcandidate.eta_2)
            Zcandidate_el_phi_2.append(Zcandidate.phi_2)
                
        out_data[self.Zel_prefix + "mass"] = Zcandidate_el_mass
        out_data[self.Zel_prefix + "pt"] = Zcandidate_el_pt
        out_data[self.Zel_prefix + "eta"] = Zcandidate_el_eta 
        out_data[self.Zel_prefix + "phi"] = Zcandidate_el_phi
        out_data[self.Zel_prefix + "dR"] = Zcandidate_el_dR
        out_data[self.Zel_prefix + "deta"] = Zcandidate_el_deta
        out_data[self.Zel_prefix + "dphi"] = Zcandidate_el_dphi
        out_data[self.Zel_prefix + "lep1_pt"] =  Zcandidate_el_lep1_pt
        out_data[self.Zel_prefix + "lep1_eta"] =  Zcandidate_el_lep1_eta
        out_data[self.Zel_prefix + "lep1_phi"] =  Zcandidate_el_lep1_phi
        out_data[self.Zel_prefix + "lep2_pt"] =  Zcandidate_el_lep2_pt
        out_data[self.Zel_prefix + "lep2_eta"] =  Zcandidate_el_lep2_eta
        out_data[self.Zel_prefix + "lep2_phi"] =  Zcandidate_el_lep2_phi
        out_data[self.Zel_prefix + "mass_2"] = Zcandidate_el_mass_2
        out_data[self.Zel_prefix + "pt_2"] = Zcandidate_el_pt_2
        out_data[self.Zel_prefix + "eta_2"] = Zcandidate_el_eta_2 
        out_data[self.Zel_prefix + "phi_2"] = Zcandidate_el_phi_2

        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True