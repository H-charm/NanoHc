import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


era_dict = {"2015": '2016preVFP_UL', "2016": '2016postVFP_UL', "2017": '2017_UL', "2018": '2018_UL'}

class ElectronSFProducer(Module, object):

    def __init__(self, year, dataset_type, **kwargs):
        self.year = year
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        correction_file = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{self.era}/electron.json.gz'
        self.corr = correctionlib.CorrectionSet.from_file(correction_file)['UL-Electron-ID-SF']

    def get_sf(self, sf_type, lep):
        if abs(lep.pdgId) != 11:
            raise RuntimeError('Input lepton is not a electron')

        wp = None
        if sf_type == 'Reco':
            wp = 'RecoBelow20' if lep.pt < 20 else 'RecoAbove20'
        elif sf_type == 'ID':
            wp = lep._wp_ID

        scale_factor = self.corr.evaluate(self.era.replace('_UL', ''), "sf", wp, lep.etaSC, lep.pt)
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
            key = 'NUM_TightID_DEN_genTracks'
        elif sf_type == 'Iso':
            key = f'NUM_{lep._wp_Iso}_DEN_TightIDandIPCut'
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
            wgtID *= self.get_sf('ID', lep)
            wgtIso *= self.get_sf('Iso', lep)

        eventWgt = wgtID * wgtIso
        self.out.fillBranch('muEffWeight', eventWgt)

        return True