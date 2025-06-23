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
            self.out.branch('elTriggerWeight', "F", limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('elHLTWeightUp', "F", limitedPrecision=10)
                self.out.branch('elHLTWeightDown', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightUp', "F", limitedPrecision=10)
                self.out.branch('elTriggerWeightDown', "F", limitedPrecision=10)

def analyze(self, event):
    """process event, return True (go to next module) or False (fail, go to next event)"""

    if not self.isMC:
        return True

    wgt = wgtUp = wgtDown = 1.0

    eff1 = eff2 = 1.0
    eff1Up = eff2Up = eff1Down = eff2Down = 1.0

    electrons = [ele for ele in event.selectedElectrons if abs(ele.pdgId) == 11]

    if len(electrons) == 2:
        el1, el2 = electrons[:2]

        if self.doSysVar:
            sf1, sf1Up, sf1Down = self.get_sf(el1)
            sf2, sf2Up, sf2Down = self.get_sf(el2)

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
            sf1 = self.get_sf(el1)
            sf2 = self.get_sf(el2)

            wgt *= sf1 * sf2
            eff1 = sf1
            eff2 = sf2

    combined_weight = eff1 + eff2 - eff1 * eff2
    self.out.fillBranch('elHLTWeight', wgt)
    self.out.fillBranch('elTriggerWeight', combined_weight)

    if self.doSysVar:
        combined_weight_up = eff1Up + eff2Up - eff1Up * eff2Up
        combined_weight_down = eff1Down + eff2Down - eff1Down * eff2Down

        self.out.fillBranch('elHLTWeightUp', wgtUp)
        self.out.fillBranch('elHLTWeightDown', wgtDown)
        self.out.fillBranch('elTriggerWeightUp', combined_weight_up)
        self.out.fillBranch('elTriggerWeightDown', combined_weight_down)

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
            self.out.branch('muTriggerWeight', 'F', limitedPrecision=10)
            if self.doSysVar:
                self.out.branch('muHLTWeightUp', 'F', limitedPrecision=10)
                self.out.branch('muHLTWeightDown', 'F', limitedPrecision=10)
                self.out.branch('muTriggerWeightUp', 'F', limitedPrecision=10)
                self.out.branch('muTriggerWeightDown', 'F', limitedPrecision=10)
            def analyze(self, event):
    """process event, return True (go to next module) or False (fail, go to next event)"""

    if not self.isMC:
        return True

    # Default weights (for multiplicative weight)
    wgt = wgtUp = wgtDown = 1.0

    # Default trigger SFs (for combined efficiency)
    eff1 = eff2 = 1.0
    eff1Up = eff2Up = eff1Down = eff2Down = 1.0

    muons = [mu for mu in event.selectedMuons if abs(mu.pdgId) == 13 and mu.pt > 26]

    if len(muons) == 2:
        mu1, mu2 = muons[:2]  # Take leading 2 muons

        if self.doSysVar:
            sf1, sf1Up, sf1Down = self.get_sf(mu1)
            sf2, sf2Up, sf2Down = self.get_sf(mu2)

            eff1 = sf1
            eff2 = sf2
            eff1Up = sf1Up
            eff2Up = sf2Up
            eff1Down = sf1Down
            eff2Down = sf2Down

            # Multiplicative event weight
            wgt *= sf1 * sf2
            wgtUp *= sf1Up * sf2Up
            wgtDown *= sf1Down * sf2Down

        else:
            sf1 = self.get_sf(mu1)
            sf2 = self.get_sf(mu2)
            eff1 = sf1
            eff2 = sf2

            wgt *= sf1 * sf2

    else:
        print('SOMETHING IS WRONG -- DID NOT FOUND 2 LEPTONS')
        eff1 = eff2 = 1.0
        if self.doSysVar:
            eff1Up = eff2Up = eff1Down = eff2Down = 1.0


    combined_weight = eff1 + eff2 - eff1 * eff2
    self.out.fillBranch('muHLTWeight', wgt)
    self.out.fillBranch('muTriggerWeight', combined_weight)

    if self.doSysVar:
        combined_weight_up = eff1Up + eff2Up - eff1Up * eff2Up
        combined_weight_down = eff1Down + eff2Down - eff1Down * eff2Down

        self.out.fillBranch('muHLTWeightUp', wgtUp)
        self.out.fillBranch('muHLTWeightDown', wgtDown)
        self.out.fillBranch('muTriggerWeightUp', combined_weight_up)
        self.out.fillBranch('muTriggerWeightDown', combined_weight_down)

    return True