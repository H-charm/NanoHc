import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}
key_dict = {"2022":     '2022Re-recoBCD',
            "2022EE":   '2022Re-recoE+PromptFG',
            "2023":     '2023PromptC',
            "2023BPix": '2023PromptD'
            }

class ElectronSFProducer(Module, object):

    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        #self.path=f'{self.year}Re-recoBCD'
        # correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electron.json.gz'
        correction_file = f'/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/ElectronSF/{self.year}/electron.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-ID-SF']

    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')

        wp = None
        if sf_type == 'Reco':
            wp = 'RecoBelow20' if lep.pt < 20 else 'Reco20to75' if 20 < lep.pt < 75 else 'RecoAbove75'
        elif sf_type == 'ID':
            wp = lep._wp_ID

        if self.year == "2022" or self.year == "2022EE":
            scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt)
        elif self.year == "2023" or self.year == "2023BPix":
            scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt, lep.phi)
        else:
            raise ValueError(f"ElectronSFProducer: Era {self.year} not supported")

        return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('elEffWeight', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtReco = 1
        wgtID = 1

        for lep in event.selectedElectrons:
            if abs(lep.pdgId) != 11:
                continue
            wgtReco *= self.get_sf('Reco', lep)
            wgtID *= self.get_sf('ID', lep)

        eventWgt = wgtReco * wgtID
        self.out.fillBranch('elEffWeight', eventWgt)

        return True


class MuonSFProducer(Module, object):

    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        correction_file_muon_Z = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_Z.json.gz' # pt binning starts from 15 
        correction_file_muon_JPsi = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_JPsi.json.gz'  # pt binning starts from 2 

        self.corr_muon_Z = correctionlib.CorrectionSet.from_file(correction_file_muon_Z)
        self.corr_muon_JPsi = correctionlib.CorrectionSet.from_file(correction_file_muon_JPsi)

    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 13:
            raise RuntimeError('Input lepton is not a muon')
        if sf_type == 'ID':
            assert lep._wp_ID == 'TightID'
            key = 'NUM_TightID_DEN_TrackerMuons'
        elif sf_type == 'Iso':
            key = f'NUM_{lep._wp_Iso}_DEN_TightID'
        if lep.pt > 15:
            scale_factor = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "nominal")
        elif lep.pt <= 15:
            scale_factor = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "nominal")
        return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('muEffWeight', "F", limitedPrecision=10)
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtID = 1
        wgtIso = 1

        for lep in event.selectedMuons:
            if abs(lep.pdgId) != 13:
                continue
            if lep.pt>15:
                wgtID *= self.get_sf('ID', lep)
                wgtIso *= self.get_sf('Iso', lep)
            else:
                wgtID *= self.get_sf('ID', lep)


        eventWgt = wgtID * wgtIso
        self.out.fillBranch('muEffWeight', eventWgt)

        return True

class MuonScaleProducer(Module, object):
    
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        self.isMC = dataset_type == "mc"
        
        if self.year == "2022":
            self.fname = "2022_schemaV2.json.gz"
        elif self.year == "2022EE":
            self.fname = "2022EE_schemaV2.json.gz"
        elif self.year == "2023":
            self.fname = "2023_schemaV2.json.gz"
        elif self.year == "2023BPix":
            self.fname = "2023BPix_schemaV2.json.gz"
        else:
            raise ValueError(f"MuonScaleProducer: Era {self.year} not supported")
        
        json_path = f"../../data/MuonScale/{self.fname}"
        
        from correctionlib import CorrectionSet
        self.corr = CorrectionSet.from_file(json_path)["cb_params"]
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('muonScale', "F", limitedPrecision=10)
    
    def analyze(self, event):
        if not self.isMC:
            return True
        
        for lep in event.selectedMuons:
            if abs(lep.pdgId) != 13:
                continue
            scale_corr = self.corr.evaluate(abs(lep.eta), float(lep.nTrackerLayers), 1)
            self.out.fillBranch('muonScale', scale_corr)
        
        return True


class ElectronScaleProducer(Module, object):
    
    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        self.isMC = dataset_type == "mc"
        
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electronSS.json.gz'
        
        from correctionlib import CorrectionSet
        self.corr = CorrectionSet.from_file(correction_file)["Scale"]
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.isMC:
            self.out = wrappedOutputTree
            self.out.branch('electronScale', "F", limitedPrecision=10)
    
    def analyze(self, event):
        if not self.isMC:
            return True
        
        for lep in event.selectedElectrons:
            if abs(lep.pdgId) != 11:
                continue
            scale_corr = self.corr.evaluate("total_correction", lep.seedGain, float(event.run), lep.eta, lep.r9, lep.pt)
            self.out.fillBranch('electronScale', scale_corr)
        
        return True
