import sys
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from ..helpers.utils import sumP4

lumi_dict = {"2022": 7.9804, "2022EE": 26.6717, "2023": 17.794, "2023BPix": 9.451}

class Zcandidate:
    def __init__(self, lep1, lep2):
        self.lep1 = lep1
        self.lep2 = lep2
        p4 = sumP4(self.lep1, self.lep2)
        self.pt = p4.Pt()
        self.eta = p4.Eta()
        self.phi = p4.Phi()
        self.mass = p4.M()
        self.dR = deltaR(lep1.eta, lep1.phi, lep2.eta, lep2.phi)
        self.deta = abs(lep1.eta - lep2.eta)
        # wrap Δφ to (−π, π], then abs
        dphi = (lep1.phi - lep2.phi) % (2*math.pi)
        if dphi > math.pi: dphi -= 2*math.pi
        self.dphi = abs(dphi)
        self.lep1_pt, self.lep1_eta, self.lep1_phi = lep1.pt, lep1.eta, lep1.phi
        self.lep2_pt, self.lep2_eta, self.lep2_phi = lep2.pt, lep2.eta, lep2.phi

class EventProducerEE(Module):

    def __init__(self, year, dataset_type, sample):
        self.year = year
        self.sample = sample
        self.dataset_type = dataset_type

        self.jet_vars = ["pt", "eta", "phi"]
        self.lep_vars = ["pt", "eta", "phi", "pdgId"]
        self.Z_vars = ["pt", "eta", "phi", "mass", "dR", "deta", "dphi",
                       "lep1_pt", "lep1_eta", "lep1_phi", "lep2_pt", "lep2_eta", "lep2_phi"]

        self.mu_prefix = "mu_"
        self.el_prefix = "el_"
        self.jet_prefix = "jet_"
        self.Zmu_prefix = "Zmu_"
        self.Zel_prefix = "Zel_"

    def beginJob(self): pass
    def endJob(self): pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = (self.dataset_type == "mc")
        self.out = wrappedOutputTree

        # branches
        for lep_var in self.lep_vars:
            self.out.branch(self.mu_prefix + lep_var, "F", 20, lenVar="nMu")
            self.out.branch(self.el_prefix + lep_var, "F", 20, lenVar="nEl")
        for jet_var in self.jet_vars:
            self.out.branch(self.jet_prefix + jet_var, "F", 20, lenVar="nJet")
        for Z_var in self.Z_vars:
            self.out.branch(self.Zmu_prefix + Z_var, "F", 20, lenVar="nZmu")
            self.out.branch(self.Zel_prefix + Z_var, "F", 20, lenVar="nZel")

        self.out.branch("HLT_pass", "O")
        self.out.branch("lumiwgt", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # PV
        if event.PV_npvsGood < 1:
            return False

        # HLT: EGamma SingleEle only (electron HLT only, for both data & MC)
        if not self._select_triggers(event):
            return False

        # Electrons first (pt>10, |η|<2.4, wp80iso)
        self._select_electrons(event)
        if len(event.selectedElectrons) != 2:
            return False

        # Event-level pT thresholds: lead>33, sublead>20
        selE = sorted(event.selectedElectrons, key=lambda x: x.pt, reverse=True)
        if not (selE[0].pt > 33 and selE[1].pt > 20):
            return False

        # Trigger matching to electrons (require ≥1 matched electron)
        if not self._trigger_match_electrons(event):
            return False

        # Muons next (pt>10, |η|<2.4, tight iso+ID), cleaned against selectedElectrons (ΔR>0.4)
        self._select_muons(event)
        if len(event.selectedMuons) != 0:
            return False  # muon_veto

        # Jets (pt>20, |η|<2.5, tightlepveto=6), cleaned from selected leptons (ΔR>0.4)
        self._select_jets(event)

        # Dielectrons: built from selectedElectrons with OS, ΔR>0.3, 60<m<120
        event.Zcandidates = []
        event.Zcandidates_el = []
        self._select_Zel_candidates(event)  # no extra pT cuts here
        if len(event.Zcandidates_el) != 1:
            return False  # one_dielectron

        # Fill outputs
        self._fill_event_info(event)
        return True

    # -------- Triggers (electron-only) --------
    def _select_triggers(self, event):
        passSingleEle = getattr(event, "HLT_Ele30_WPTight_Gsf", False)
        passTrigger = bool(passSingleEle)
        self.out.fillBranch("HLT_pass", passTrigger)
        return passTrigger

    # -------- Trigger matching  --------
    def _trigger_match_electrons(self, event, max_dr=0.1):
        """
        Mirrors:
          - id == 11
          - filterbit bit 1 (1<<1)
          - trigobj pt > 28
          - ΔR(e, trigobj) < 0.1
        Require >= 1 matched electron.
        """
        try:
            trigobjs = Collection(event, "TrigObj")
        except Exception:
            return False

        # select trigger objects passing Ele30_WPTight_Gsf filters
        bit = (1 << 1)
        to_sel = []
        for to in trigobjs:
            if abs(getattr(to, "id", 0)) != 11:  # electron
                continue
            if getattr(to, "pt", 0.0) <= 28.0:
                continue
            if (getattr(to, "filterBits", 0) & bit) == 0:
                continue
            to_sel.append(to)
        if not to_sel:
            return False

        # require any selected electron to match one of the selected trig objects
        for el in event.selectedElectrons:
            for to in to_sel:
                if deltaR(el.eta, el.phi, to.eta, to.phi) < max_dr:
                    return True
        return False

    # -------- Object selections --------
    def _select_electrons(self, event):
        event.selectedElectrons = []
        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt <= 10: continue
            if abs(el.eta) >= 2.4: continue
            if abs(el.dxy) >= 0.5 or abs(el.dz) >= 1.0: continue
            if not getattr(el, "mvaIso_WP80", False): continue
            el._wp_ID = 'wp80iso'
            event.selectedElectrons.append(el)

    def _select_muons(self, event):
        event.selectedMuons = []
        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt <= 10: continue
            if abs(mu.eta) >= 2.4: continue
            if abs(mu.dxy) >= 0.5 or abs(mu.dz) >= 1.0: continue
            # tight iso + tight ID
            if mu.pfRelIso03_all > 0.15: continue
            if not getattr(mu, "tightId", False): continue

            # clean vs *selected* electrons only (ΔR>0.4)
            overlaps = False
            for el in getattr(event, "selectedElectrons", []):
                if deltaR(mu.eta, mu.phi, el.eta, el.phi) <= 0.4:
                    overlaps = True
                    break
            if overlaps: continue

            mu._wp_ID = 'TightID'
            mu._wp_Iso = 'TightPFIso'
            event.selectedMuons.append(mu)

    def _select_jets(self, event):
        event.selectedJets = []
        jets = Collection(event, "Jet")
        for jet in jets:
            if jet.pt <= 20: continue
            if abs(jet.eta) >= 2.5: continue
            if jet.jetId != 6: continue  # tightlepveto (Run-3 mapping)

            # clean from selected leptons (ΔR>0.4)
            overlaps = False
            for lep in getattr(event, "selectedElectrons", []):
                if deltaR(lep.eta, lep.phi, jet.eta, jet.phi) <= 0.4:
                    overlaps = True; break
            if not overlaps:
                for lep in getattr(event, "selectedMuons", []):
                    if deltaR(lep.eta, lep.phi, jet.eta, jet.phi) <= 0.4:
                        overlaps = True; break
            if overlaps: continue

            event.selectedJets.append(jet)

    def _select_Zel_candidates(self, event):
        # All OS pairs from selectedElectrons; ΔR>0.3; 60 < m < 120
        els = getattr(event, "selectedElectrons", [])
        if len(els) < 2:
            return
        for i in range(len(els)):
            for j in range(i+1, len(els)):
                el1, el2 = els[i], els[j]
                # Opposite sign (use charge)
                if el1.charge * el2.charge >= 0:
                    continue
                # ΔR > 0.3
                if deltaR(el1.eta, el1.phi, el2.eta, el2.phi) <= 0.3:
                    continue
                Zcand_el = Zcandidate(el1, el2) if el1.pt >= el2.pt else Zcandidate(el2, el1)
                if 60.0 < Zcand_el.mass < 120.0:
                    event.Zcandidates_el.append(Zcand_el)
                    event.Zcandidates.append(Zcand_el)

    # (unchanged) branch filling
    def _fill_event_info(self, event):
        out_data = {}
        out_data["lumiwgt"] = lumi_dict[self.year]

        leptons_pt_sorted = sorted(event.selectedElectrons + event.selectedMuons,
                                   key=lambda p: p.pt, reverse=True)

        el_pt, el_eta, el_phi, el_pdgId = [], [], [], []
        mu_pt, mu_eta, mu_phi, mu_pdgId = [], [], [], []

        for lep in leptons_pt_sorted:
            if abs(lep.pdgId) == 11:
                el_pt.append(lep.pt); el_eta.append(lep.eta); el_phi.append(lep.phi); el_pdgId.append(lep.pdgId)
            elif abs(lep.pdgId) == 13:
                mu_pt.append(lep.pt); mu_eta.append(lep.eta); mu_phi.append(lep.phi); mu_pdgId.append(lep.pdgId)

        out_data[self.el_prefix + "pt"] = el_pt
        out_data[self.el_prefix + "eta"] = el_eta
        out_data[self.el_prefix + "phi"] = el_phi
        out_data[self.el_prefix + "pdgId"] = el_pdgId
        out_data[self.mu_prefix + "pt"] = mu_pt
        out_data[self.mu_prefix + "eta"] = mu_eta
        out_data[self.mu_prefix + "phi"] = mu_phi
        out_data[self.mu_prefix + "pdgId"] = mu_pdgId

        jet_pt, jet_eta, jet_phi = [], [], []
        for jet in event.selectedJets:
            jet_pt.append(jet.pt); jet_eta.append(jet.eta); jet_phi.append(jet.phi)
        out_data[self.jet_prefix + "pt"] = jet_pt
        out_data[self.jet_prefix + "eta"] = jet_eta
        out_data[self.jet_prefix + "phi"] = jet_phi

        Zm = []; Zpt = []; Zeta = []; Zphi = []
        ZdR = []; Zdeta = []; Zdphi = []
        Zl1pt = []; Zl1eta = []; Zl1phi = []
        Zl2pt = []; Zl2eta = []; Zl2phi = []
        for Z in event.Zcandidates_el:
            Zm.append(Z.mass); Zpt.append(Z.pt); Zeta.append(Z.eta); Zphi.append(Z.phi)
            ZdR.append(Z.dR); Zdeta.append(Z.deta); Zdphi.append(Z.dphi)
            Zl1pt.append(Z.lep1_pt); Zl1eta.append(Z.lep1_eta); Zl1phi.append(Z.lep1_phi)
            Zl2pt.append(Z.lep2_pt); Zl2eta.append(Z.lep2_eta); Zl2phi.append(Z.lep2_phi)

        out_data[self.Zel_prefix + "mass"] = Zm
        out_data[self.Zel_prefix + "pt"] = Zpt
        out_data[self.Zel_prefix + "eta"] = Zeta
        out_data[self.Zel_prefix + "phi"] = Zphi
        out_data[self.Zel_prefix + "dR"] = ZdR
        out_data[self.Zel_prefix + "deta"] = Zdeta
        out_data[self.Zel_prefix + "dphi"] = Zdphi
        out_data[self.Zel_prefix + "lep1_pt"] = Zl1pt
        out_data[self.Zel_prefix + "lep1_eta"] = Zl1eta
        out_data[self.Zel_prefix + "lep1_phi"] = Zl1phi
        out_data[self.Zel_prefix + "lep2_pt"] = Zl2pt
        out_data[self.Zel_prefix + "lep2_eta"] = Zl2eta
        out_data[self.Zel_prefix + "lep2_phi"] = Zl2phi

        # muon Z list kept for structure compatibility (likely empty in Zee)
        out_data[self.Zmu_prefix + "mass"] = []
        out_data[self.Zmu_prefix + "pt"] = []
        out_data[self.Zmu_prefix + "eta"] = []
        out_data[self.Zmu_prefix + "phi"] = []
        out_data[self.Zmu_prefix + "dR"] = []
        out_data[self.Zmu_prefix + "deta"] = []
        out_data[self.Zmu_prefix + "dphi"] = []
        out_data[self.Zmu_prefix + "lep1_pt"] = []
        out_data[self.Zmu_prefix + "lep1_eta"] = []
        out_data[self.Zmu_prefix + "lep1_phi"] = []
        out_data[self.Zmu_prefix + "lep2_pt"] = []
        out_data[self.Zmu_prefix + "lep2_eta"] = []
        out_data[self.Zmu_prefix + "lep2_phi"] = []

        for key, val in out_data.items():
            self.out.fillBranch(key, val)