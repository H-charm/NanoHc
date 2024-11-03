from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from ..helpers.triggerHelper import passTrigger
import ROOT
import math
import itertools
from ..helpers.utils import sumP4
ROOT.PyConfig.IgnoreCommandLineOptions = True

class Zcandidate:
    pass


class Hcandidate:
    pass

        
class ZZcandidate:
    pass


class HcTreeProducer(Module):
    
    def __init__(self, year, dataset_type):
        self.year = year
        self.dataset_type = dataset_type

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        
        self.out = wrappedOutputTree
        
        ## define lepton branches
        self.out.branch("mu_pt", "F", 20, lenVar="nMu")
        self.out.branch("mu_eta", "F", 20, lenVar="nMu")
        self.out.branch("mu_phi", "F", 20, lenVar="nMu")
        self.out.branch("el_pt", "F", 20, lenVar="nEl")
        self.out.branch("el_eta", "F", 20, lenVar="nEl")
        self.out.branch("el_phi", "F", 20, lenVar="nEl")
        self.out.branch("lep_pt", "F", 20, lenVar="nLeptons")
        self.out.branch("lep_eta", "F", 20, lenVar="nLeptons")
        self.out.branch("lep_phi", "F", 20, lenVar="nLeptons")
        
        ## define jet branches
        self.out.branch("n_jets", "I")
        self.out.branch("ak4_bdisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_cvbdisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_cvldisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_gvudsdisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_pt", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_eta", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_phi", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_mass", "F", 20, lenVar="n_ak4")
        if self.isMC: self.out.branch("ak4_hadronFlavour", "F", 20, lenVar="n_ak4")
        
        ## define trigger branches
        self.out.branch("passTriggers", "O")
        
        ## Zcandidates
        self.out.branch("Zcandidate_mass", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_pt", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_eta", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_phi", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_onshell_mass", "F", 20, lenVar="n_Zcandidates_onshell")
        self.out.branch("Zcandidate_offshell_mass", "F", 20, lenVar="n_Zcandidates_offshell")
  
        ## Hcandidates
        self.out.branch("Hcandidate_mass", "F", 20, lenVar="n_Hcandidates")
        self.out.branch("Hcandidate_pt", "F", 20, lenVar="n_Hcandidates")
        self.out.branch("Hcandidate_eta", "F", 20, lenVar="n_Hcandidates")
        self.out.branch("Hcandidate_phi", "F", 20, lenVar="n_Hcandidates")
        
        if self.isMC:
            self.out.branch("l1PreFiringWeight", "F", limitedPrecision=10)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if self._select_triggers(event) is False:
            return False

        self._select_muons(event)
        self._select_electrons(event)  
        event.selectedLeptons = event.selectedMuons + event.selectedElectrons
        if len(event.selectedLeptons) < 4: return False

        self._select_jets(event)
        if len(event.selectedJets) == 0:
            return False
        
        self._select_Z_candidates(event)
        if len(event.Zcandidates) < 2:
            return False
        
        self._select_ZZ_candidates(event)
        if len(event.ZZcandidates) > 1: return False # for now keep only events with 1 ZZ candidate (we need MELA if we have more)
        
        self._select_H_candidates(event)
        
        self._fill_event_info(event)

        return True
    
    def _select_triggers(self, event):

        out_data = {}

        if self.year == "2016APV":
            pass #FIX
        if self.year == "2016":
            pass #FIX
        if self.year == "2017":
            out_data["passTriggers"] = passTrigger(event, [
                'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL',
                'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL',
                'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 
                'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8', 
                'HLT_TripleMu_12_10_5', 
                'HLT_TripleMu_10_5_5_D2', 
                'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 
                'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
                'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ', 
                'HLT_Mu8_DiEle12_CaloIdL_TrackIdL', 
                'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ', 
                'HLT_Ele35_WPTight_Gsf_v', 
                'HLT_Ele38_WPTight_Gsf_v', 
                'HLT_Ele40_WPTight_Gsf_v', 
                'HLT_IsoMu27',
            ])
        elif self.year == "2018":
            out_data["passTriggers"] = passTrigger(event, [
                'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'HLT_DoubleEle25_CaloIdL_MW',
                'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8',
                'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 
                'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
                'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
                'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ', 
                'HLT_Ele32_WPTight_Gsf', 
                'HLT_IsoMu24', 
            ])
        else:
            print(f"Year {self.year} not found")
            
        # apply trigger selections on data
        if not self.isMC: 
            if not out_data['passTriggers']:
                return False

        for key in out_data:
            self.out.fillBranch(key, out_data[key])

        return True

    def _select_Z_candidates(self, event):

        event.Zcandidates = []
        
        lepton_pairs = list(itertools.combinations(event.selectedLeptons, 2))
        
        for lepton_pair in lepton_pairs:
            
            # we need same flavor and opposite charge
            if (lepton_pair[0].pdgId + lepton_pair[1].pdgId) != 0:
                continue
            
            # let's put negative charged lepton always as lep1 (useful later)
            lep1, lep2 = (lepton_pair[0], lepton_pair[1]) if lepton_pair[0].pdgId < 0 else (lepton_pair[1], lepton_pair[0])
            
            Vboson = sumP4(lep1, lep2)
            if Vboson.M() < 12 or Vboson.M() > 120:
                continue
                        
            Zcand = Zcandidate()
            Zcand.pt = Vboson.Pt()
            Zcand.mass = Vboson.M()
            Zcand.eta = Vboson.Eta()
            Zcand.phi = Vboson.Phi()
            Zcand.lep1 = lep1
            Zcand.lep2 = lep2
        
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
            
            if Z1.mass < 40: continue  
            if Z1.lep1.pt < 20: continue
            if Z1.lep2.pt < 10: continue
            
            leptons = [Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2]
            lepton_pairs = list(itertools.combinations(leptons, 2))
            
            dr_ll_values = [] # between each of the four leptons 
            m_ll_values = [] # between each of the four leptons (opposite-sign pairs regardless of flavors)                     
            for lepton_pair in lepton_pairs:
                lep1 = lepton_pair[0]
                lep2 = lepton_pair[1]
                
                dr = math.sqrt( (lep1.eta - lep2.eta)**2 + (lep1.phi - lep2.phi)**2 ) 
                dr_ll_values.append(dr)
                
                if lep1.pdgId * lep2.pdgId < 0:
                    m_ll = sumP4(lep1, lep2).M()
                    m_ll_values.append(m_ll)
                    
            if any(dr < 0.02 for dr in dr_ll_values): continue
            if any(m_ll < 4 for m_ll in m_ll_values): continue
            
            ## check alternative pairing (4e or 4mu)
            if abs(Z1.lep1.pdgId) == abs(Z1.lep2.pdgId) == abs(Z2.lep1.pdgId) == abs(Z2.lep2.pdgId):

                Za_lep1 = Z2.lep1
                Za_lep2 = Z1.lep2
                
                Zb_lep1 = Z1.lep1
                Zb_lep2 = Z2.lep2
                
                mZ = 91.1876
                
                ## reject if |mZa - mZ| < |mZ1 - mZ| and mZb < 12
                Za_mass = sumP4(Za_lep1, Za_lep2).M()
                Zb_mass = sumP4(Zb_lep1, Zb_lep2).M()
                if ( abs(Za_mass - mZ) < abs(Z1.mass - mZ) ) and Zb_mass < 12: continue

                ## reject if inv mass of 4-lepton system < 70              
                m_4l = sumP4(Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2).M()
                if m_4l < 70: continue
                
            ZZcand = ZZcandidate()
            ZZcand.Z1 = Z1
            ZZcand.Z2 = Z2
            event.ZZcandidates.append(ZZcand)      
            
    def _select_H_candidates(self, event):
        
        event.Hcandidates = []
                    
        for ZZcand in event.ZZcandidates:
            
            Z1 = ZZcand.Z1
            Z2 = ZZcand.Z2
            Vboson = sumP4(Z1, Z2)
            
            Hcand = Hcandidate()
            Hcand.pt = Vboson.Pt()
            Hcand.mass = Vboson.M()
            Hcand.eta = Vboson.Eta()
            Hcand.phi = Vboson.Phi()
            Hcand.Z1 = Z1
            Hcand.Z2 = Z2
                        
            event.Hcandidates.append(Hcand)
            
    def _select_muons(self, event):

        event.selectedMuons = []

        muons = Collection(event, "Muon")
        
        for mu in muons:
            
            if mu.pt > 5 and abs(mu.eta) < 2.4 and mu.dxy < 0.5 and mu.dz < 1 and abs(mu.sip3d) < 4 and mu.pfRelIso03_all < 0.35 and mu.isPFcand:
                mu._wp_ID = 'TightID'
                mu._wp_Iso = 'LooseRelIso'
                event.selectedMuons.append(mu)
                          
    def _select_electrons(self, event):

        event.selectedElectrons = []

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.5 and el.dxy < 0.5 and el.dz < 1 and abs(el.sip3d) < 4:
                el._wp_ID = 'wp90iso'
                event.selectedElectrons.append(el)

    def _select_jets(self, event):

        event.selectedJets = []

        jets = Collection(event, "Jet")
        photons = Collection(event, "Photon")

        for jet in jets:
            if jet.pt <= 15 and abs(jet.eta) >= 2.5:
                continue
            
            jet_isolated = True
            for lep in event.selectedLeptons:
                dR_jet_lep = math.sqrt( (lep.eta - jet.eta)**2 + (lep.phi - jet.phi)**2 ) 
                if dR_jet_lep <= 0.4:
                    jet_isolated = False
            if not jet_isolated:
                continue
                    
            photon_isolated = True
            for photon in photons:
                dR_jet_photon = math.sqrt( (photon.eta - jet.eta)**2 + (photon.phi - jet.phi)**2 ) 
                if dR_jet_photon <= 0.4:
                    photon_isolated = False
            if not photon_isolated:
                continue                   
            
            event.selectedJets.append(jet)
        
    def _fill_event_info(self, event):
        out_data = {}
        
        ## leptons
        leptons_pt_sorted = sorted(event.selectedLeptons, key=lambda particle: particle.pt, reverse=True)

        lep_pt = []
        lep_eta = []
        lep_phi = []
        el_pt = []
        el_eta = []
        el_phi = []
        mu_pt = []
        mu_eta = []
        mu_phi = []
        for lep in leptons_pt_sorted:
            lep_pt.append(lep.pt)
            lep_eta.append(lep.eta)
            lep_phi.append(lep.phi)
            if abs(lep.pdgId) == 11: # el
                el_pt.append(lep.pt)
                el_eta.append(lep.eta)
                el_phi.append(lep.phi)
            if abs(lep.pdgId) == 13: # mu
                mu_pt.append(lep.pt)
                mu_eta.append(lep.eta)
                mu_phi.append(lep.phi)            
            
        out_data["lep_pt"] = lep_pt
        out_data["lep_eta"] = lep_eta
        out_data["lep_phi"] = lep_phi
        out_data["el_pt"] = el_pt
        out_data["el_eta"] = el_eta
        out_data["el_phi"] = el_phi       
        out_data["mu_pt"] = mu_pt
        out_data["mu_eta"] = mu_eta
        out_data["mu_phi"] = mu_phi
                    
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
            ak4_bdisc.append(jet.btagDeepFlavB)
            ak4_cvbdisc.append(jet.btagDeepFlavCvB)
            ak4_cvldisc.append(jet.btagDeepFlavCvL)
            ak4_gvudsdisc.append(jet.btagDeepFlavQG)
            ak4_pt.append(jet.pt)
            ak4_eta.append(jet.eta)
            ak4_phi.append(jet.phi)
            ak4_mass.append(jet.mass)
            if self.isMC: ak4_hadronFlavour.append(jet.hadronFlavour)
        
        out_data["ak4_bdisc"] = ak4_bdisc
        out_data["ak4_cvbdisc"] = ak4_cvbdisc
        out_data["ak4_cvldisc"] = ak4_cvldisc 
        out_data["ak4_gvudsdisc"] = ak4_gvudsdisc 
        out_data["ak4_pt"] = ak4_pt 
        out_data["ak4_eta"] = ak4_eta 
        out_data["ak4_phi"] = ak4_phi 
        out_data["ak4_mass"] = ak4_mass 
        if self.isMC: out_data["ak4_hadronFlavour"] = ak4_hadronFlavour 
           
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
                
        out_data["Zcandidate_mass"] = Zcandidate_mass
        out_data["Zcandidate_pt"] = Zcandidate_pt
        out_data["Zcandidate_eta"] = Zcandidate_eta 
        out_data["Zcandidate_phi"] = Zcandidate_phi 
        out_data["Zcandidate_onshell_mass"] = Zcandidate_onshell_mass
        out_data["Zcandidate_offshell_mass"] = Zcandidate_offshell_mass

        ## H candidates
        Hcandidate_mass = []
        Hcandidate_pt = []
        Hcandidate_eta = []
        Hcandidate_phi = []
        for Hcandidate in event.Hcandidates:
            Hcandidate_mass.append(Hcandidate.mass)
            Hcandidate_pt.append(Hcandidate.pt)
            Hcandidate_eta.append(Hcandidate.eta)
            Hcandidate_phi.append(Hcandidate.phi)
            
        out_data["Hcandidate_mass"] = Hcandidate_mass
        out_data["Hcandidate_pt"] = Hcandidate_pt
        out_data["Hcandidate_eta"] = Hcandidate_eta 
        out_data["Hcandidate_phi"] = Hcandidate_phi 
                
        if self.isMC:
            out_data["l1PreFiringWeight"] = event.L1PreFiringWeight_Nom                
                
        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True