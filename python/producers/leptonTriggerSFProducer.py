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

class ElectronHLTSFProducer(Module, object):

    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]
        #self.path=f'{self.year}Re-recoBCD'
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electronHlt.json.gz'
        # correction_file = f'../../data/ElectronSF/{self.year}/electron.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-HLT-SF']

    def get_sf(self, lep):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')

        #WPs: "HLT_SF_Ele30_LooseID","HLT_SF_Ele30_MediumID","HLT_SF_Ele30_TightID","HLT_SF_Ele30_MVAiso80ID", "HLT_SF_Ele30_MVAiso90ID"
        wp = 'HLT_SF_Ele30_MVAiso90ID'

        scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt)
            if self.doSysVar:
                scale_factor_up = self.corr.evaluate(key_dict[self.year], "sfup", wp, lep.etaSC, lep.pt)
                scale_factor_down = self.corr.evaluate(key_dict[self.year], "sfdown", wp, lep.etaSC, lep.pt)

        return scale_factor, scale_factor_up, scale_factor_down if self.doSysVar else scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('elHLTWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elHLTWeightUp', "F", limitedPrecision=10)
                self.out.branch('elHLTWeightDown', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgt = wgtUp = wgtDown = 1

        for lep in event.selectedElectrons:
            if abs(lep.pdgId) != 11:
                continue

            if self.doSysVar:
                sf, sfUp, sfDown = self.get_sf(lep)

                wgt *= sf
                wgtUp *= sfUp
                wgtDown *= sfDown
            else:
                wgt *= self.get_sf(lep)

        self.out.fillBranch('elHLTWeight', wgt)

        if self.doSysVar:
            self.out.fillBranch('elHLTWeightUp', wgtUp)
            self.out.fillBranch('elHLTWeightDown', wgtDown)

        return True

# 'nominal', 'stat', 'syst', 'systup', and 'systdown' (systup = syst + stat, systdown = syst - stat)
class MuonHLTSFProducer(Module, object):

    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]
        correction_file_muon_Z = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_Z.json.gz' # pt binning starts from 15 
        #correction_file_muon_JPsi = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{self.era}/muon_JPsi.json.gz'  # pt binning starts from 2 

        self.corr_muon_Z = correctionlib.CorrectionSet.from_file(correction_file_muon_Z)
        #self.corr_muon_JPsi = correctionlib.CorrectionSet.from_file(correction_file_muon_JPsi)

    def get_sf(self, lep):
        if abs(lep.pdgId) != 13:
            raise RuntimeError('Input lepton is not a muon')
        key = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
        if lep.pt > 26:
            scale_factor = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "nominal")
            if self.doSysVar:
                scale_factor_up = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "systup")
                scale_factor_down = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "systdown")
        # elif lep.pt <= 15:
        #     scale_factor = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "nominal")
        #     if self.doSysVar:
        #         scale_factor_up = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "systup")
        #         scale_factor_down = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "systdown")
        return scale_factor, scale_factor_up, scale_factor_down if self.doSysVar else scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('muHLTWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muHLTWeightUp', 'F', limitedPrecision=10)
                self.out.branch('muHLTWeightDown', 'F', limitedPrecision=10)
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgt = wgtUp = wgtDown = 1

        for lep in event.selectedMuons:
            if abs(lep.pdgId) != 13:
                continue
            if lep.pt>26:
                if self.doSysVar:
                    sf, sfUp, sfDown = self.get_sf(lep)
                    wgt *= sf
                    wgtUp *= sfUp
                    wgtDown *= sfDow

                else:
                    wgt *= self.get_sf(lep)


        
        self.out.fillBranch('muHLTWeight', wgt)

        if self.doSysVar:
            self.out.fillBranch('muHLTWeightUp', wgtUp)
            self.out.fillBranch('muHLTWeightDown', wgtDown)
        return True