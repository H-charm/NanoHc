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

class ElectronTRGSFProducer(Module, object):

    def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.doSysVar = doSysVar
        self.era = era_dict[self.year]
        #self.path=f'{self.year}Re-recoBCD'
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electronHlt.json.gz'
        # correction_file = f'../../data/ElectronSF/{self.year}/electron.json.gz'
        self.corr_mc = correctionlib.CorrectionSet.from_file(correction_file)['Electron-HLT-McEff']
        self.corr_data = correctionlib.CorrectionSet.from_file(correction_file)['Electron-HLT-DataEff']
        
    def get_sf(self, lep, corr_type):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')
        if corr_type == 'mc':
            corr = self.corr_mc
        else:
            corr = self.corr_data

        #WPs: "HLT_SF_Ele30_LooseID","HLT_SF_Ele30_MediumID","HLT_SF_Ele30_TightID","HLT_SF_Ele30_MVAiso80ID", "HLT_SF_Ele30_MVAiso90ID"
        wp = 'HLT_SF_Ele30_MVAiso80ID'
        scale_factor_up = scale_factor_down = scale_factor = 1
        if lep.pt > 25:
            scale_factor = corr.evaluate(key_dict[self.year], "nom", wp, lep.etaSC, lep.pt)
            if self.doSysVar:
                scale_factor_up = corr.evaluate(key_dict[self.year], "up", wp, lep.etaSC, lep.pt)
                scale_factor_down = corr.evaluate(key_dict[self.year], "down", wp, lep.etaSC, lep.pt)

        if self.doSysVar:
            return scale_factor, scale_factor_up, scale_factor_down
        else:
            return scale_factor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = True if self.dataset_type == "mc" else False
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('elHLTWeight', "F", limitedPrecision=10)
            self.out.branch('elTriggerWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elHLTWeightUp', "F", limitedPrecision=10)
                self.out.branch('elHLTWeightDown', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightUp', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightDown', "F", limitedPrecision=10)

            self.out.branch('elHLTWeightData', "F", limitedPrecision=10)
            self.out.branch('elTriggerWeightData', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elHLTWeightUpData', "F", limitedPrecision=10)
                self.out.branch('elHLTWeightDownData', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightUpData', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightDownData', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if not self.isMC:
            return True

        corr_type = ['mc', 'data']
        for corr in corr_type:

            wgt = wgtUp = wgtDown = 1.0

            eff1 = eff2 = 1.0
            eff1Up = eff2Up = eff1Down = eff2Down = 1.0

            if len(event.selectedElectrons) == 2:
                el1, el2 = event.selectedElectrons[:2]

                if self.doSysVar:
                    sf1, sf1Up, sf1Down = self.get_sf(el1, corr)
                    sf2, sf2Up, sf2Down = self.get_sf(el2, corr)

                    wgt *= sf1 * sf2
                    wgtUp *= sf1Up * sf2Up
                    wgtDown *= sf1Down * sf2Down

                    eff1 = sf1
                    eff2 = sf2
                    eff1Up = sf1Up
                    eff2Up = sf2Up
                    eff1Down = sf1Down
                    eff2Down = sf2Down
                else:
                    sf1 = self.get_sf(el1, 'mc')
                    sf2 = self.get_sf(el2, 'mc')

                    wgt *= sf1 * sf2
                    eff1 = sf1
                    eff2 = sf2
            else:
                return True 

            combined_weight = eff1 + eff2 - eff1 * eff2

            if corr == 'mc':
                self.out.fillBranch('elHLTWeight', wgt)
                self.out.fillBranch('elTriggerWeight', combined_weight)

                if self.doSysVar:
                    combined_weight_up = eff1Up + eff2Up - eff1Up * eff2Up
                    combined_weight_down = eff1Down + eff2Down - eff1Down * eff2Down

                    self.out.fillBranch('elHLTWeightUp', wgtUp)
                    self.out.fillBranch('elHLTWeightDown', wgtDown)
                    self.out.fillBranch('elTriggerWeightUp', combined_weight_up)
                    self.out.fillBranch('elTriggerWeightDown', combined_weight_down)
            else:
                self.out.fillBranch('elHLTWeightData', wgt)
                self.out.fillBranch('elTriggerWeightData', combined_weight)

                if self.doSysVar:
                    combined_weight_up = eff1Up + eff2Up - eff1Up * eff2Up
                    combined_weight_down = eff1Down + eff2Down - eff1Down * eff2Down

                    self.out.fillBranch('elHLTWeightUpData', wgtUp)
                    self.out.fillBranch('elHLTWeightDownData', wgtDown)
                    self.out.fillBranch('elTriggerWeightUpData', combined_weight_up)
                    self.out.fillBranch('elTriggerWeightDownData', combined_weight_down)

        return True

# class ElectronHLTSFProducer(Module, object):

#     def __init__(self, year, dataset_type, doSysVar=False, **kwargs):
#         self.year = year
#         self.dataset_type = dataset_type
#         self.doSysVar = doSysVar
#         self.era = era_dict[self.year]
#         #self.path=f'{self.year}Re-recoBCD'
#         correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electronHlt.json.gz'
#         # correction_file = f'../../data/ElectronSF/{self.year}/electron.json.gz'
#         self.corr = correctionlib.CorrectionSet.from_file(correction_file)['Electron-HLT-SF']

#     def get_sf(self, lep):
#         if abs(lep.pdgId) != 11:
#             raise RuntimeError('Input lepton is not a electron')

#         #WPs: "HLT_SF_Ele30_LooseID","HLT_SF_Ele30_MediumID","HLT_SF_Ele30_TightID","HLT_SF_Ele30_MVAiso80ID", "HLT_SF_Ele30_MVAiso90ID"
#         wp = 'HLT_SF_Ele30_MVAiso80ID'
#         scale_factor_up = scale_factor_down = scale_factor = 1
#         if lep.pt > 25:
#             scale_factor = self.corr.evaluate(key_dict[self.year], "sf", wp, lep.etaSC, lep.pt)
#             if self.doSysVar:
#                 scale_factor_up = self.corr.evaluate(key_dict[self.year], "sfup", wp, lep.etaSC, lep.pt)
#                 scale_factor_down = self.corr.evaluate(key_dict[self.year], "sfdown", wp, lep.etaSC, lep.pt)

#         if self.doSysVar:
#             return scale_factor, scale_factor_up, scale_factor_down
#         else:
#             return scale_factor

#     def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#         self.isMC = True if self.dataset_type == "mc" else False
#         if self.isMC:
#             self.out = wrappedOutputTree

#             self.out.branch('elHLTSFWeight', "F", limitedPrecision=10)
#             self.out.branch('elTriggerSFWeight', "F", limitedPrecision=10)
#             if self.doSysVar:
#                 self.out.branch('elHLTSFWeightUp', "F", limitedPrecision=10)
#                 self.out.branch('elHLTSFWeightDown', "F", limitedPrecision=10)
#                 self.out.branch('elTriggerSFWeightUp', "F", limitedPrecision=10)
#                 self.out.branch('elTriggerSFWeightDown', "F", limitedPrecision=10)


#     def analyze(self, event):
#         """process event, return True (go to next module) or False (fail, go to next event)"""
#         if not self.isMC:
#             return True

#         wgt = wgtUp = wgtDown = 1.0

#         eff1 = eff2 = 1.0
#         eff1Up = eff2Up = eff1Down = eff2Down = 1.0

#         if len(event.selectedElectrons) == 2:
#             el1, el2 = event.selectedElectrons[:2]

#             if self.doSysVar:
#                 sf1, sf1Up, sf1Down = self.get_sf(el1)
#                 sf2, sf2Up, sf2Down = self.get_sf(el2)

#                 wgt *= sf1 * sf2
#                 wgtUp *= sf1Up * sf2Up
#                 wgtDown *= sf1Down * sf2Down

#                 eff1 = sf1
#                 eff2 = sf2
#                 eff1Up = sf1Up
#                 eff2Up = sf2Up
#                 eff1Down = sf1Down
#                 eff2Down = sf2Down
#             else:
#                 sf1 = self.get_sf(el1)
#                 sf2 = self.get_sf(el2)

#                 wgt *= sf1 * sf2
#                 eff1 = sf1
#                 eff2 = sf2
#         else:
#             return True 

#         combined_weight = eff1 + eff2 - eff1 * eff2

#         self.out.fillBranch('elHLTSFWeight', wgt)
#         self.out.fillBranch('elTriggerSFWeight', combined_weight)

#         if self.doSysVar:
#             combined_weight_up = eff1Up + eff2Up - eff1Up * eff2Up
#             combined_weight_down = eff1Down + eff2Down - eff1Down * eff2Down

#             self.out.fillBranch('elHLTSFWeightUp', wgtUp)
#             self.out.fillBranch('elHLTSFWeightDown', wgtDown)
#             self.out.fillBranch('elTriggerSFWeightUp', combined_weight_up)
#             self.out.fillBranch('elTriggerSFWeightDown', combined_weight_down)

#         return True