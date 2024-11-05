#########################################################################
# Module to add some additional lepton variables as a branch to nanoAOD #
#########################################################################

# imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
# from pathlib import Path

# import nanoAODTools
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class LeptonVariablesModule(Module):

    def __init__( self, variables=['all'] ):
        ### intializer
        # input arguments:
        # - variables: list of variables names to add
        #   (default: add all variables defined here)
        self.variables = variables
        self.defined_variables = ['jetPtRatio', 'jetBTagDeepFlavor']
        if 'all' in self.variables:
            self.variables = self.defined_variables
        # check provided variables
        for variable in self.variables:
            if variable not in self.defined_variables:
                msg = 'ERROR in LeptonVariablesModule:'
                msg += ' variable {} not recognized; options are {}'.format(
                       variable, self.defined_variables)
                raise Exception(msg)
        print('Initialized a LeptonVariablesModule with following parameters:')
        print('  - variables:')
        for var in self.variables: print('    - {}'.format(var))

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for variable in self.variables:
            self.out.branch('Electron_{}'.format(variable), "F", lenVar='nElectron')
            self.out.branch('Muon_{}'.format(variable), "F", lenVar='nMuon')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        ### process a single event
        # (always return True as this module performs no selection)

        # get electrons and muons
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")

        # jetPtRatio
        if 'jetPtRatio' in self.variables:
            electron_jetptratio = [1./(el.jetRelIso+1) for el in electrons]
            muon_jetptratio = [1./(mu.jetRelIso+1) for mu in muons]
            self.out.fillBranch('Electron_jetPtRatio', electron_jetptratio)
            self.out.fillBranch('Muon_jetPtRatio', muon_jetptratio)

        # jetBTagDeepFlavor
        if 'jetBTagDeepFlavor' in self.variables:
            electron_jdf = ([ (jets[el.jetIdx].btagDeepFlavB if el.jetIdx>=0 else 0)
                              for el in electrons ])
            muon_jdf = ([ (jets[mu.jetIdx].btagDeepFlavB if mu.jetIdx>=0 else 0)
                          for mu in muons ])
            self.out.fillBranch('Electron_jetBTagDeepFlavor', electron_jdf)
            self.out.fillBranch('Muon_jetBTagDeepFlavor', muon_jdf)

        return True