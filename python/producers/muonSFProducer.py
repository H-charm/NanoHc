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

# 'nominal', 'stat', 'syst', 'systup', and 'systdown' (systup = syst + stat, systdown = syst - stat)
class MuonSFProducer(Module, object):

    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
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
        scale_factor_up = scale_factor_down = 1
        if lep.pt > 15:
            scale_factor = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "nominal")
            if self.doSysVar:
                scale_factor_up = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "systup")
                scale_factor_down = self.corr_muon_Z[key].evaluate(abs(lep.eta), lep.pt, "systdown")
        elif lep.pt <= 15:
            scale_factor = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "nominal")
            if self.doSysVar:
                scale_factor_up = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "systup")
                scale_factor_down = self.corr_muon_JPsi[key].evaluate(abs(lep.eta), lep.pt, "systdown")
        if self.doSysVar:
            return scale_factor, scale_factor_up, scale_factor_down
        else:
            return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('muEffWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muEffWeightUp', 'F', limitedPrecision=10)
                self.out.branch('muEffWeightDown', 'F', limitedPrecision=10)
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtID = wgtIDUp = wgtIDDown = 1
        wgtIso = wgtIsoUp = wgtIsoDown = 1

        for lep in event.selectedMuons:
            if abs(lep.pdgId) != 13:
                continue
            if self.doSysVar:
                if lep.pt>15:
                    sfID, sfIDUp, sfIDDown = self.get_sf('ID', lep)
                    sfIso, sfIsoUp, sfIsoDown = self.get_sf('Iso', lep)

                    wgtID *= sfID
                    wgtIDUp *= sfIDUp
                    wgtIDDown *= sfIDDown
                    wgtIso *= sfIso
                    wgtIsoUp *= sfIsoUp
                    wgtIsoDown *= sfIsoDown

                else:
                    sfID, sfIDUp, sfIDDown = self.get_sf('ID', lep)
                    wgtID *= sfID
                    wgtIDUp *= sfIDUp
                    wgtIDDown *= sfIDDown
            else:
                if lep.pt>15:
                    sfID = self.get_sf('ID', lep)
                    sfIso = self.get_sf('Iso', lep)
                    wgtID *= sfID
                    wgtIso *= sfIso
                else:
                    sfID = self.get_sf('ID', lep)
                    wgtID *= sfID

        eventWgt = wgtID * wgtIso
        self.out.fillBranch('muEffWeight', eventWgt)

        if self.doSysVar:
            eventWgtUp = wgtIsoUp * wgtIDUp
            eventWgtDown = wgtIsoDown * wgtIDDown
            self.out.fillBranch('muEffWeightUp', eventWgtUp)
            self.out.fillBranch('muEffWeightDown', eventWgtDown)
        return True