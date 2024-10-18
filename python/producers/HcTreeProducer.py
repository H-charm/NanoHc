from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from ..helpers.triggerHelper import passTrigger
import ROOT
import math
import itertools
from ..helpers.utils import sumP4
ROOT.PyConfig.IgnoreCommandLineOptions = True

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
        self.out.branch("n_muons", "I")
        self.out.branch("n_electrons", "I")
        for lep_idx in range(4):
            for lep_var in ["pt","eta","phi","pdgId","charge"]:
                self.out.branch("lep" + str(lep_idx+1) + "_" + lep_var, "F")
        
        ## define jet branches
        self.out.branch("n_jets", "I")
        self.out.branch("ak4_bdisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_cvbdisc", "F", 20, lenVar="n_ak4")
        self.out.branch("ak4_cvldisc", "F", 20, lenVar="n_ak4")
        
        ## define trigger branches
        self.out.branch("passTriggers", "O")
        
        ## Zcandidates
        self.out.branch("Zcandidate_mass", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_pt", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_eta", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_phi", "F", 20, lenVar="n_Zcandidates")
        self.out.branch("Zcandidate_onshell_mass", "F")
        self.out.branch("Zcandidate_offshell_mass", "F")
  
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

        ## trigger selection
        if self._selectTriggers(event) is False:
            return False

        ## lepton selection
        self._selectMuons(event)
        self._selectElectrons(event)  
        
        if len(event.Muons) + len(event.Electrons) != 4:
            return False
        
        if len(event.Muons) not in [0,2,4] or len(event.Electrons) not in [0,2,4]: # we can only have 0,2,4 muons / electrons
            return False
        
        event.Leptons = event.Muons + event.Electrons
        leptons_charge_sum = 0
        for lepton in event.Leptons:
            leptons_charge_sum += lepton.charge
        if leptons_charge_sum != 0:
            return False
        
        ## select jets
        if self._selectJets(event) is False:
            return False      
        
        if self._selectZcandidates(event) is False:
            return False
        
        self._find_onshell_and_offshell_Zcandidates(event)   
        
        self._selectHcandidates(event)
        
        self._fillEventInfo(event)

        return True
    
    def _selectTriggers(self, event):

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

    def _selectZcandidates(self, event):

        event.Zcandidates = []
        
        leptons = event.Muons + event.Electrons
        lepton_pairs = list(itertools.combinations(leptons, 2))
        
        for lepton_pair in lepton_pairs:
            lep1 = lepton_pair[0]
            lep2 = lepton_pair[1]
            
            if lep1.pdgId * lep2.pdgId > 0: # we need OS 
                continue
            if (lep1.pdgId + lep2.pdgId) != 0: # we need same flavor
                continue
            
            dR_lep1_lep2 = math.sqrt( (lep1.eta - lep2.eta)**2 + (lep1.phi - lep2.phi)**2 )
            if dR_lep1_lep2 <= 0.02:
                continue
            
            Vboson = sumP4(lep1, lep2)
            if Vboson.M() < 10 or Vboson.M() > 120:
                continue
            
            event.Zcandidates.append(Vboson)
            
        if len(event.Zcandidates) != 2: # we need 2 Z candidates for ZZ candidate
            return False
        else:
            return True
        
    
    def _find_onshell_and_offshell_Zcandidates(self, event):
        
        Z1 = event.Zcandidates[0]
        Z2 = event.Zcandidates[1]
        
        mZ = 91.1876
        
        ## minimal |mZcandidate - mZ| is considered on-shell
        d_mZ1_mZ = abs( Z1.M() - mZ )
        d_mZ2_mZ = abs( Z2.M() - mZ )
        
        onshell_idx = 0 if d_mZ1_mZ < d_mZ2_mZ else 1
        offshell_idx = 1 if d_mZ1_mZ < d_mZ2_mZ else 0
        
        out_data = {}
        
        out_data["Zcandidate_onshell_mass"] = event.Zcandidates[onshell_idx].M()
        out_data["Zcandidate_offshell_mass"] = event.Zcandidates[offshell_idx].M()
        
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
                    
        
    def _selectHcandidates(self, event):
        
        event.Hcandidates = []
        
        # these are P4 vectors
        Z1 = event.Zcandidates[0]
        Z2 = event.Zcandidates[1]
        
        Hcandidate = Z1 + Z2
        
        event.Hcandidates.append(Hcandidate)
        
        
    def _selectMuons(self, event):

        event.Muons = []

        muons = Collection(event, "Muon")
        
        for mu in muons:
            
            if mu.pt > 30 and abs(mu.eta) < 2.4 and mu.dxy < 0.5 and mu.dz < 1 and abs(mu.sip3d) < 4 and mu.pfRelIso03_all < 0.35 and mu.tightId:
                mu._wp_ID = 'TightID'
                mu._wp_Iso = 'LooseRelIso'
                event.Muons.append(mu)
                
                
    def _selectElectrons(self, event):

        event.Electrons = []

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.5 and el.dxy < 0.5 and el.dz < 1 and abs(el.sip3d) < 4:
                el._wp_ID = 'wp90iso'
                event.Electrons.append(el)

    def _selectJets(self, event):

        event.Jets = []

        jets = Collection(event, "Jet")
        for jet in jets:
            if jet.pt <= 15 and abs(jet.eta) >= 2.5:
                continue
            
            jet_isolated = True
            for lep in event.Leptons:
                dR_jet_lep = math.sqrt( (lep.eta - jet.eta)**2 + (lep.phi - jet.phi)**2 ) 
                if dR_jet_lep <= 0.4:
                    jet_isolated = False
            if not jet_isolated:
                continue
                    
            event.Jets.append(jet)
        
        if len(event.Jets) == 0:
            return False
        

    def _fillEventInfo(self, event):
        out_data = {}
        
        ## leptons
        out_data["n_muons"] = len(event.Muons)
        out_data["n_electrons"] = len(event.Electrons)

        pt_sorted_leptons = sorted(event.Leptons, key=lambda obj: obj.pt, reverse=True)
        for lep_idx in range(4):
            for lep_var in ["pt","phi","eta","pdgId","charge"]:
                out_data["lep" + str(lep_idx+1) + "_" + lep_var] = getattr(pt_sorted_leptons[lep_idx],lep_var)
        
        
        ## jets  
        out_data["n_jets"] = len(event.Jets)
        ak4_bdisc = []
        ak4_cvbdisc = []
        ak4_cvldisc = []
        
        for jet in event.Jets:
            ak4_bdisc.append(jet.btagDeepFlavB)
            ak4_cvbdisc.append(jet.btagDeepFlavCvB)
            ak4_cvldisc.append(jet.btagDeepFlavCvL)
        
        out_data["ak4_bdisc"] = ak4_bdisc
        out_data["ak4_cvbdisc"] = ak4_cvbdisc
        out_data["ak4_cvldisc"] = ak4_cvldisc 
           
        ## Z candidates
        Zcandidate_mass = []
        Zcandidate_pt = []
        Zcandidate_eta = []
        Zcandidate_phi = []
        for Zcandidate in event.Zcandidates:
            Zcandidate_mass.append(Zcandidate.M())
            Zcandidate_pt.append(Zcandidate.Pt())
            Zcandidate_eta.append(Zcandidate.Eta())
            Zcandidate_phi.append(Zcandidate.Phi())
            
        out_data["Zcandidate_mass"] = Zcandidate_mass
        out_data["Zcandidate_pt"] = Zcandidate_pt
        out_data["Zcandidate_eta"] = Zcandidate_eta 
        out_data["Zcandidate_phi"] = Zcandidate_phi 

        ## H candidates
        Hcandidate_mass = []
        Hcandidate_pt = []
        Hcandidate_eta = []
        Hcandidate_phi = []
        for Hcandidate in event.Hcandidates:
            Hcandidate_mass.append(Hcandidate.M())
            Hcandidate_pt.append(Hcandidate.Pt())
            Hcandidate_eta.append(Hcandidate.Eta())
            Hcandidate_phi.append(Hcandidate.Phi())
            
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