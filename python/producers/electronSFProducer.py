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

    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]
        #self.path=f'{self.year}Re-recoBCD'
        # correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electron.json.gz'
        correction_file = f'../../data/ElectronSF/{self.year}/electron.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-ID-SF']

    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')

        wp = None
        scale_factor_up = scale_factor_down = 1
        if sf_type == 'Reco':
            wp = 'RecoBelow20' if lep.pt < 20 else 'Reco20to75' if 20 < lep.pt < 75 else 'RecoAbove75'
        elif sf_type == 'ID':
            wp = lep._wp_ID

        if self.year == "2022" or self.year == "2022EE":
            scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt)
            if self.doSysVar:
                scale_factor_up = self.corr.evaluate(key_dict[self.year], "sfup", wp, lep.etaSC, lep.pt)
                scale_factor_down = self.corr.evaluate(key_dict[self.year], "sfdown", wp, lep.etaSC, lep.pt)
        elif self.year == "2023" or self.year == "2023BPix":
            scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt, lep.phi)
            if self.doSysVar:
                scale_factor_up = self.corr.evaluate(key_dict[self.year], "sfup", wp, lep.etaSC, lep.pt, lep.phi)
                scale_factor_down = self.corr.evaluate(key_dict[self.year], "sfdown", wp, lep.etaSC, lep.pt, lep.phi)
        else:
            raise ValueError(f"ElectronSFProducer: Era {self.year} not supported")

        if self.doSysVar:
            return scale_factor, scale_factor_up, scale_factor_down
        else:
            return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('elEffWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elEffWeightUp', "F", limitedPrecision=10)
                self.out.branch('elEffWeightDown', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

        wgtReco = wgtRecoUp = wgtRecoDown = 1
        wgtID = wgtIDUp = wgtIDDown = 1

        for lep in event.selectedElectrons:
            if abs(lep.pdgId) != 11:
                continue

            if self.doSysVar:
                sfReco, sfRecoUp, sfRecoDown = self.get_sf('Reco', lep)
                sfID, sfIDUp, sfIDDown = self.get_sf('ID', lep)

                wgtReco *= sfReco
                wgtRecoUp *= sfRecoUp
                wgtRecoDown *= sfRecoDown

                wgtID *= sfID
                wgtIDUp *= sfIDUp
                wgtIDDown *= sfIDDown
            else:
                sfReco = self.get_sf('Reco', lep)
                sfID = self.get_sf('ID', lep)
                wgtReco *= sfReco
                wgtID *= sfID

        eventWgt = wgtReco * wgtID
        self.out.fillBranch('elEffWeight', eventWgt)

        if self.doSysVar:
            eventWgtUp = wgtRecoUp * wgtIDUp
            eventWgtDown = wgtRecoDown * wgtIDDown
            self.out.fillBranch('elEffWeightUp', eventWgtUp)
            self.out.fillBranch('elEffWeightDown', eventWgtDown)

        return True

