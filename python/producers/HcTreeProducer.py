from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class HcTreeProducer(Module):
    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("n_muons", "I")
        self.out.branch("n_electrons", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self._selectMuons(event)
        self._selectElectrons(event)
        
        if len(event.Muons) + len(event.Electrons) != 4:
            return False
        if len(event.Muons) not in [0,2,4] or len(event.Electrons) not in [0,2,4]: # we can only have 0,2,4 muons / electrons
            return False

        self._fillEventInfo(event)

        return True

    def _selectMuons(self, event):

        event.Muons = []

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 5 and abs(mu.eta) < 2.4 and mu.dxy < 0.5 and mu.dz < 1 and abs(mu.sip3d) < 4 and mu.pfRelIso03_all < 0.35 and mu.tightId:
                event.Muons.append(mu)

    def _selectElectrons(self, event):

        event.Electrons = []

        electrons = Collection(event, "Electron")
        for el in electrons:
            if el.pt > 7 and abs(el.eta) < 2.5 and el.dxy < 0.5 and el.dz < 1 and abs(el.sip3d) < 4:
                event.Electrons.append(el)

    def _selectJets(self, event):

        event.Jets = []

        jets = Collection(event, "Jet")
        for jet in jets:
            if jet.pt > 15 and abs(jet.eta) < 2.5:
                event.Jets.append(jet)

    def _fillEventInfo(self, event):
        out_data = {}
        
        out_data["n_muons"] = len(event.Muons)
        out_data["n_electrons"] = len(event.Electrons)
        
        for key in out_data:
            self.out.fillBranch(key, out_data[key])
            
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
HcTreeProducerModule = lambda: HcTreeProducer()