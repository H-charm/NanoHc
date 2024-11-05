##################################################################
# Module to add the TOP lepton MVA scores as a branch to nanoAOD #
##################################################################
# based on this implementation:
# https://github.com/HephyAnalysisSW/Analysis/blob/UL/Tools/python/mvaTOPreader.py
# note: this module requires the variables jetPtRatio and jetBTagDeepFlavor
#       which are not in the standard nanoAOD but are set by LeptonVariablesModule;
#       hence always run TopLeptonMvaModule after LeptonVariablesModule!

# imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
# from pathlib import Path
import numpy as np
import xgboost as xgb

# import nanoAODTools
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class TopLeptonMvaModule(Module):

    def __init__( self, year, version,
                  verbose=True, variablename='mvaTOP' ):
        ### intializer
        # input arguments:
        # - year: data taking year, used to retrieve the correct MVA weights
        # - version: MVA version, choose from 'ULv1' or 'ULv2'
        # - verbose: print more or less output
        # - variablename: name of the variable to add as a branch,
        #   e.g. the default 'mvaTOP' will add branches
        #   Electron_mvaTOP and Muon_mvaTOP
        # - crab: whether this is run in a crab job or locally
        #   (which impacts where to find the BDT weight files)
        self.year = year
        self.version = version
        self.electronvarname = 'Electron_{}'.format(variablename)
        self.muonvarname = 'Muon_{}'.format(variablename)

        # check arguments
        if year not in ['2016PreVFP','2016PostVFP','2017','2018']:
            msg = 'ERROR in TopLeptonMvaModule:'
            msg += ' year {} not recognized.'.format(year)
            raise Exception(msg)
        if version not in ['ULv1', 'ULv2']:
            msg = 'ERROR in TopLeptonMvaModule:'
            msg += ' year {} not recognized.'.format(year)
            raise Exception(msg)

        # set directory and file names
        weightdir = os.path.join(os.path.dirname(__file__), '../../data/leptonmva/weights')
        if not os.path.exists(weightdir):
            # for CRAB submission, the data directory is copied to the working directory
            weightdir = 'data/leptonmva/weights'
        if not os.path.exists(weightdir):
            raise Exception('ERROR: weight directory not found.')
        weightfile = 'TOP'
        if self.version == 'ULv2': weightfile += 'v2'
        diryear = year.replace('20','')
        if year=='2016PreVFP': diryear = '16APV'
        if year=='2016PostVFP': diryear = '16'
        weightfile += 'UL' + diryear + '_XGB.weights.bin'
        elweightfile = os.path.join(weightdir, 'el_' + weightfile)
        muweightfile = os.path.join(weightdir, 'mu_' + weightfile)

        # do printouts if requested
        if verbose:
            print('Initializing a TopLeptonMvaModule with following properties:')
            print('  - year: {}'.format(year))
            print('  - version: {}'.format(version))
            print('  - electron weights file: {}'.format(elweightfile))
            print('  - muon weights file: {}'.format(muweightfile))

        # check if weight files exist
        for f in [elweightfile, muweightfile]:
            if not os.path.exists(f):
                msg = 'ERROR in TopLeptonMvaModule:'
                msg += ' file {} does not exist.'.format(f)
                raise Exception(msg)

        # load weights
        self.electronmva = xgb.Booster()
        self.electronmva.load_model(elweightfile)
        self.muonmva = xgb.Booster()
        self.muonmva.load_model(muweightfile)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.electronvarname, "F", lenVar='nElectron')
        self.out.branch(self.muonvarname, "F", lenVar='nMuon')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getElectronMvaScore(self, electron):
        ### get the MVA score for an electron
        features = ([[
            electron.pt, 
            electron.eta, 
            electron.jetNDauCharged,
            electron.miniPFRelIso_chg,
            electron.miniPFRelIso_all - electron.miniPFRelIso_chg,
            electron.jetPtRelv2,
            electron.jetPtRatio,
            electron.pfRelIso03_all,
            electron.jetBTagDeepFlavor,
	    electron.sip3d,
	    np.log(abs(electron.dxy)),
	    np.log(abs(electron.dz)),
            electron.mvaFall17V2noIso,
        ]])
        if self.version=='ULv2': features[0].append( electron.lostHits )
        features = np.array(features)
        fmatrix = xgb.DMatrix(features, nthread=1)
        score = self.electronmva.predict(fmatrix)[0]
        return score

    def getMuonMvaScore(self, muon):
        ### get the MVA score for an muon
        features = ([[
            muon.pt,
            muon.eta,
            muon.jetNDauCharged,
            muon.miniPFRelIso_chg,
            muon.miniPFRelIso_all - muon.miniPFRelIso_chg,
            muon.jetPtRelv2,
            muon.jetPtRatio,
            muon.pfRelIso03_all,
            muon.jetBTagDeepFlavor,
            muon.sip3d,
            np.log(abs(muon.dxy)),
            np.log(abs(muon.dz)),
            muon.segmentComp
        ]])
        features = np.array(features)
        fmatrix = xgb.DMatrix(features, nthread=1)
        score = self.muonmva.predict(fmatrix)[0]
        return score

    def analyze(self, event):
        ### process a single event
        # (always return True as this module performs no selection)

        # get electrons and muons
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")

        # calculate the mva score
        electron_scores = [self.getElectronMvaScore(el) for el in electrons]
        muon_scores = [self.getMuonMvaScore(mu) for mu in muons]
 
        # fill branches
        self.out.fillBranch(self.electronvarname, electron_scores)
        self.out.fillBranch(self.muonvarname, muon_scores)

        return True