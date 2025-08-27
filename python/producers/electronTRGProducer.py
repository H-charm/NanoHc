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

class ElectronTriggerProducer(Module, object):

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

            self.out.branch('elTriggerWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elTriggerWeightUp', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightDown', "F", limitedPrecision=10)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if not self.isMC:
            return True

        eff1 = eff2 = 1.0
        eff1Up = eff2Up = eff1Down = eff2Down = 1.0

        corr_type = ['mc', 'data']

        combined_weight = {
            "mc": {"nom": 1.0, "up": 1.0, "down": 1.0},
            "data": {"nom": 1.0, "up": 1.0, "down": 1.0},
        }

        if len(event.selectedElectrons) == 2:
            el1, el2 = event.selectedElectrons[:2]
            for corr in corr_type:

                if self.doSysVar:
                    sf1, sf1Up, sf1Down = self.get_sf(el1, corr)
                    sf2, sf2Up, sf2Down = self.get_sf(el2, corr)

                    eff1, eff2 = sf1, sf2
                    eff1Up, eff2Up = sf1Up, sf2Up
                    eff1Down, eff2Down = sf1Down, sf2Down
                else:
                    sf1 = self.get_sf(el1, 'mc')
                    sf2 = self.get_sf(el2, 'mc')

                    eff1, eff2 = sf1, sf2
                
                combined_weight[corr]["nom"] = eff1 + eff2 - eff1 * eff2
                if self.doSysVar:
                    combined_weight[corr]["up"] = eff1Up + eff2Up - eff1Up * eff2Up
                    combined_weight[corr]["down"] = eff1Down + eff2Down - eff1Down * eff2Down
        else:
            return True 

        if combined_weight["mc"]["nom"] != 0:
            sf = combined_weight["data"]["nom"]/ combined_weight["mc"]["nom"]
        else:
            sf = 1.0

        self.out.fillBranch('elTriggerWeight', sf)
        
        if self.doSysVar:
            sf_up = (
            combined_weight["data"]["up"] / combined_weight["mc"]["up"] if combined_weight["mc"]["up"] != 0 else 1.0
        )
            sf_down = (
            combined_weight["data"]["down"] / combined_weight["mc"]["down"] if combined_weight["mc"]["down"] != 0 else 1.0
        )
            self.out.fillBranch("elTriggerWeightUp", sf_up)
            self.out.fillBranch("elTriggerWeightDown", sf_down)

        return True