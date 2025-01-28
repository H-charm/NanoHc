from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from ..helpers.triggerHelper import passTrigger
import ROOT
import math
import itertools
from ..helpers.utils import sumP4
ROOT.PyConfig.IgnoreCommandLineOptions = True

class Zcandidate:
    
    def __init__(self,lep1,lep2):
        self.lep1 = lep1
        self.lep2 = lep2
        self.pt = sumP4(self.lep1, self.lep2).Pt()
        self.eta = sumP4(self.lep1, self.lep2).Eta()
        self.phi = sumP4(self.lep1, self.lep2).Phi()
        self.mass = sumP4(self.lep1, self.lep2).M()


class ZZcandidate:

    def __init__(self,Z1,Z2):
        self.Z1 = Z1
        self.Z2 = Z2
        self.pt = sumP4(self.Z1, self.Z2).Pt()
        self.eta = sumP4(self.Z1, self.Z2).Eta()
        self.phi = sumP4(self.Z1, self.Z2).Phi()
        self.mass = sumP4(self.Z1, self.Z2).M()


class BaselineProducer(Module):
    
    def __init__(self, year, dataset_type):
        self.year = year
        self.dataset_type = dataset_type
        self.lep_vars = ["pt","eta","phi"]
        self.jet_vars = ["pt","eta","phi","mass","bdisc","cvbdisc","cvldisc","gvudsdisc"]
        self.jet_vars_mc = ["hadronFlavour"]
        self.Z_vars = ["pt","eta","phi","mass","onshell_mass","offshell_mass"]
        self.H_vars = ["pt","eta","phi","mass"]        
        self.mu_prefix = "mu_"
        self.el_prefix = "el_"
        self.lep_prefix = "lep_"
        self.jet_prefix = "jet_"
        self.Z_prefix = "Z_"
        self.ZZ_prefix = "ZZ_"
        self.H_prefix = "H_"

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        
        self.out = wrappedOutputTree
        
        ## define lepton branches
        for lep_var in self.lep_vars:
            self.out.branch(self.mu_prefix + lep_var, "F", 20, lenVar="nMu")
            self.out.branch(self.el_prefix + lep_var, "F", 20, lenVar="nEl")
            self.out.branch(self.lep_prefix + lep_var, "F", 20, lenVar="nLep")
        
        ## define jet branches
        for jet_var in self.jet_vars:
            self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar="nJet")
        if self.isMC: 
            for jet_var in self.jet_vars_mc:
                self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar="nJet")
        
        ## define trigger branches
        self.out.branch("passTriggers", "O")
        
        ## Zcandidates
        for Z_var in self.Z_vars:
            self.out.branch(self.Z_prefix + Z_var, "F", 20, lenVar="nZ")
  
        ## Hcandidates
        for H_var in self.H_vars:
            self.out.branch(self.H_prefix + H_var, "F", 20, lenVar="nH")
                    
        # if self.isMC:
        #     self.out.branch("l1PreFiringWeight", "F", limitedPrecision=10)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if event.PV_npvsGood < 1: return False

        # apply trigger selections on data
        if not self.isMC: 
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
        # if len(event.ZZcandidates) > 1: return False # for now keep only events with 1 ZZ candidate (we need MELA if we have more)
        
        self._select_H_candidates(event)
        
        self._fill_event_info(event)

        return True
        
    def _select_triggers(self, event):

        out_data = {}

        if self.year == "2022":
            out_data["passTriggers"] = passTrigger(event, [
                'HLT_Ele30_WPTight_Gsf',
                'HLT_IsoMu24',
                'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'HLT_DoubleEle25_CaloIdL_MW',
                'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8',
                'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ',
                'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ',
                'HLT_TripleMu_10_5_5_DZ',
                'HLT_TripleMu_12_10_5',
            ])

        # if self.year == "2016APV":
        #     pass #FIX
        # if self.year == "2016":
        #     pass #FIX
        # if self.year == "2017":
        #     out_data["passTriggers"] = passTrigger(event, [
        #         'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
        #         'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL',
        #         'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL',
        #         'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 
        #         'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8', 
        #         'HLT_TripleMu_12_10_5', 
        #         'HLT_TripleMu_10_5_5_D2', 
        #         'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 
        #         'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
        #         'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
        #         'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
        #         'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ', 
        #         'HLT_Mu8_DiEle12_CaloIdL_TrackIdL', 
        #         'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ', 
        #         'HLT_Ele35_WPTight_Gsf_v', 
        #         'HLT_Ele38_WPTight_Gsf_v', 
        #         'HLT_Ele40_WPTight_Gsf_v', 
        #         'HLT_IsoMu27',
        #     ])
        # elif self.year == "2018":
        #     out_data["passTriggers"] = passTrigger(event, [
        #         'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
        #         'HLT_DoubleEle25_CaloIdL_MW',
        #         'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8',
        #         'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 
        #         'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
        #         'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 
        #         'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ', 
        #         'HLT_Ele32_WPTight_Gsf', 
        #         'HLT_IsoMu24', 
        #     ])
        else:
            print(f"Year {self.year} not found")
                    
        if not out_data['passTriggers']:
            return False

        for key in out_data:
            self.out.fillBranch(key, out_data[key])

        return True

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
                        
            Zcand = Zcandidate(lep1,lep2)
            
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
   
            if Z1.mass < 40: continue  
            
            leptons = [Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2] 
            
            ## two DISTINCT leptons must pass pt > 10 and pt > 20
            at_least_one_passed_pt20 = False
            at_least_one_passed_pt10 = False
            for lep in leptons:
                if lep.pt > 20: at_least_one_passed_pt20 = True
                elif lep.pt > 10: at_least_one_passed_pt10 = True
            if not (at_least_one_passed_pt10 and at_least_one_passed_pt20): continue
            
            lepton_pairs = list(itertools.combinations(leptons, 2)) # 6 combinations
            
            dr_ll_values = [] # between each of the four leptons (Ghost removal: all pairs)
            m_ll_values = [] # between each of the four leptons (QCD suppression: opposite-sign pairs and same flavour)                     
            for lepton_pair in lepton_pairs:
                lep1 = lepton_pair[0]
                lep2 = lepton_pair[1]
                            
                dr = math.sqrt( (lep1.eta - lep2.eta)**2 + (lep1.phi - lep2.phi)**2 ) 
                dr_ll_values.append(dr)
                
                if (lep1.pdgId + lep2.pdgId) == 0:
                    m_ll = sumP4(lep1, lep2).M()
                    m_ll_values.append(m_ll)
                    
            if any(dr < 0.02 for dr in dr_ll_values): continue
            if any(m_ll < 4 for m_ll in m_ll_values): continue
            
            ## smart cut: check alternative pairing (4e or 4mu)
            if abs(Z1.lep1.pdgId) == abs(Z1.lep2.pdgId) == abs(Z2.lep1.pdgId) == abs(Z2.lep2.pdgId):

                ## define Za as the one closest to Z mass, and Zb as the other pair
                Ztemp1 = Zcandidate(Z1.lep1, Z2.lep2)
                Ztemp2 = Zcandidate(Z2.lep1, Z1.lep2)

                mZ = 91.1876
                Za, Zb = (Ztemp1, Ztemp2) if ( abs(Ztemp1.mass - mZ) < abs(Ztemp2.mass - mZ) ) else (Ztemp2, Ztemp1)
                                
                ## reject if |mZa - mZ| < |mZ1 - mZ| and mZb < 12
                if ( abs(Za.mass - mZ) < abs(Z1.mass - mZ) ) and Zb.mass < 12: continue

            ## reject if inv mass of 4-lepton system < 105             
            m_4l = sumP4(Z1.lep1, Z1.lep2, Z2.lep1, Z2.lep2).M()
            if m_4l < 105: continue

            ## reject if inv mass of 4-lepton system > 140
            if m_4l > 140: continue
                
            ZZcand = ZZcandidate(Z1,Z2)
            event.ZZcandidates.append(ZZcand)      
            
    def _select_H_candidates(self, event):
        
        event.Hcandidates = []
        
        ## for now keep all ZZ canidates as Higgs candidates            
        for ZZcand in event.ZZcandidates:
            
            Hcand = ZZcand
                        
            event.Hcandidates.append(Hcand)
        
    ## taken from here https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/nanoZZ4lAnalysis.py    
    def _select_muons(self, event):

        event.selectedMuons = []

        muons = Collection(event, "Muon")
        
        for mu in muons:
            
            passMuID = mu.isPFcand or (mu.highPtId>0 and mu.pt>200)
            
            if mu.pt > 5 and abs(mu.eta) < 2.4 and mu.dxy < 0.5 and mu.dz < 1 and abs(mu.sip3d) < 4 and mu.pfRelIso03_all < 0.35 and passMuID and (mu.isGlobal or (mu.isTracker and mu.nStations>0)):
                mu._wp_ID = 'TightID'
                mu._wp_Iso = 'LoosePFIso'
                event.selectedMuons.append(mu)
            
    ## taken from here https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/nanoZZ4lAnalysis.py       
    def _select_electrons(self, event):

        event.selectedElectrons = []

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.5 and el.dxy < 0.5 and el.dz < 1 and abs(el.sip3d) < 4:
                el._wp_ID = 'wp90iso'
                
                ## https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/getEleBDTCut.py#L22-L31
                if abs(el.etaSC) < 0.8:
                    if el.pt < 10:
                        #if el.mvaIso < 0.9128577458: continue
                        if el.mvaIso < 1.6339: continue
                    else:
                        #if el.mvaIso < 0.1559788054: continue
                        if el.mvaIso < 0.3685: continue
                elif 0.8 < abs(el.etaSC) < 1.479:
                    if el.pt < 10:
                        #if el.mvaIso < 0.9056792368: continue
                        if el.mvaIso < 1.5499: continue
                    else:
                        #if el.mvaIso < 0.0273863727: continue
                        if el.mvaIso < 0.2662: continue                    
                else: # |el.etaSC| > 1.479
                    if el.pt < 10:
                        #if el.mvaIso < 0.9439440575: continue
                        if el.mvaIso < 2.0629: continue
                    else:
                        #if el.mvaIso < -0.5532483665: continue
                        if el.mvaIso < -0.5444: continue                         
                                    
                event.selectedElectrons.append(el)

    def _select_jets(self, event):

        event.selectedJets = []

        jets = Collection(event, "Jet")
        FsrPhotons = Collection(event, "FsrPhoton")

        for jet in jets:
            if jet.pt <= 15 and abs(jet.eta) >= 2.5:
                continue
            if abs(jet.phi) > math.pi: # Introduced due to jetvetomaps corrections
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
            
        out_data[self.lep_prefix + "pt"] = lep_pt
        out_data[self.lep_prefix + "eta"] = lep_eta
        out_data[self.lep_prefix + "phi"] = lep_phi
        out_data[self.el_prefix + "pt"] = el_pt
        out_data[self.el_prefix + "eta"] = el_eta
        out_data[self.el_prefix + "phi"] = el_phi       
        out_data[self.mu_prefix + "pt"] = mu_pt
        out_data[self.mu_prefix + "eta"] = mu_eta
        out_data[self.mu_prefix + "phi"] = mu_phi
                    
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
        
        out_data[self.jet_prefix + "bdisc"] = ak4_bdisc
        out_data[self.jet_prefix + "cvbdisc"] = ak4_cvbdisc
        out_data[self.jet_prefix + "cvldisc"] = ak4_cvldisc 
        out_data[self.jet_prefix + "gvudsdisc"] = ak4_gvudsdisc 
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
            
        out_data[self.H_prefix + "mass"] = Hcandidate_mass
        out_data[self.H_prefix + "pt"] = Hcandidate_pt
        out_data[self.H_prefix + "eta"] = Hcandidate_eta 
        out_data[self.H_prefix + "phi"] = Hcandidate_phi 
                
        # if self.isMC:
        #     out_data["l1PreFiringWeight"] = event.L1PreFiringWeight_Nom                
                
        ## fill all branches
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True