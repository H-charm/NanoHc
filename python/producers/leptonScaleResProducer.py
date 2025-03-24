###
# Add branches for electron scale and resolution corrections.
# 
# Example: 2022EE MC smearing+uncertainties:
# eleSS = eleScaleRes("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS.json.gz", "Scale", "Smearing", True)
# 
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
from math import pi
import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ..helpers.utils import polarP4
import numpy as np
import random
from .roccor import RoccoR


class eleScaleRes(Module):
    def __init__(self, year, dataset_type, overwritePt=False):
        self.year = year
        self.dataset_type = dataset_type
        self.overwritePt = overwritePt
        self.isMC = self.dataset_type == "mc"

        file_map = {
            2016: 'Legacy2016_07Aug2017_FineEtaR9_v3_ele_unc',
            2017: 'Run2017_17Nov2017_v1_ele_unc',
            2018: 'Run2018_Step2Closure_CoarseEtaR9Gain_v2',
        }
        
        if self.year in file_map:
            correction_file = f'EgammaAnalysis/ElectronTools/data/ScalesSmearings/{file_map[self.year]}'
            print(f'Loading egamma scale and smearing file from {correction_file}')
            self.eleCorr = ROOT.EnergyScaleCorrection(correction_file, ROOT.EnergyScaleCorrection.ECALELF)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt:
            self.out.branch("Electron_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_uncorrected_pt", "F", lenVar="nElectron")
        else:
            self.out.branch("Electron_corrected_pt", "F", lenVar="nElectron")
        if self.isMC:
            self.out.branch("Electron_scaleUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_scaleDn_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearDn_pt", "F", lenVar="nElectron")

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        pt_corr, pt_smear_up, pt_smear_dn, pt_scale_up, pt_scale_dn = [], [], [], [], []
        
        for ele in electrons:
            if self.isMC:
                self.rnd.SetSeed(self.rndSeed(event, ele))
                smearNrSigma = self.rnd.Gaus(0, 1)

                run = event.run
                abseta = abs(ele.eta + ele.deltaEtaSC)
                et = ele.energy / math.cosh(abseta)
                r9 = ele.r9
                seed = ele.seedGain

                smearCorr = self.eleCorr.getSmearCorr(run, et, abseta, r9, seed)
                if smearCorr:
                    smear = smearCorr.sigma(et)
                    smearUp = smearCorr.sigma(et, 1.0, 0.0)
                    smearDn = smearCorr.sigma(et, -1.0, 0.0)
                else:
                    smear, smearUp, smearDn = 0, 0, 0

                corrUpDelta = abs((smearUp - smear) * smearNrSigma)
                corrDnDelta = abs((smearDn - smear) * smearNrSigma)
                smearErr = max(corrUpDelta, corrDnDelta)
                scaleErr = self.eleCorr.scaleCorrUncert(run, et, abseta, r9, seed, ROOT.std.bitset(
                    ROOT.EnergyScaleCorrection.kErrNrBits)(ROOT.EnergyScaleCorrection.kErrStatSystGain))
                err = math.hypot(smearErr, scaleErr)

                pt_corr.append(ele.pt * (1 + smear * smearNrSigma))
                pt_smear_up.append(ele.pt * (1 + smear + err))
                pt_smear_dn.append(ele.pt * (1 + smear - err))
                pt_scale_up.append(ele.pt * (1 + scaleErr))
                pt_scale_dn.append(ele.pt * (1 - scaleErr))
            else:
                scale = self.eleCorr.scaleCorrUncert(run, et, abseta, r9, seed, ROOT.std.bitset(
                    ROOT.EnergyScaleCorrection.kErrNrBits)(ROOT.EnergyScaleCorrection.kErrStatSystGain))
                pt_corr.append(scale * ele.pt)

        if self.overwritePt:
            self.out.fillBranch("Electron_uncorrected_pt", [ele.pt for ele in electrons])
            self.out.fillBranch("Electron_pt", pt_corr)
        else:
            self.out.fillBranch("Electron_corrected_pt", pt_corr)

        if self.isMC:
            self.out.fillBranch("Electron_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Electron_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Electron_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Electron_scaleDn_pt", pt_scale_dn)

        return True

def rndSeed(event, obj):
    """
    Generates a deterministic random seed for reproducibility.
    """
    return (event.run << 20) + (event.luminosityBlock << 10) + event.event + int(obj.eta / 0.01)

def mk_safe(fct, *args):
    try:
        return fct(*args)
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in arg
               for arg in e.args):
            print(
                'WARNING: catching exception and returning -1. Exception arguments: %s'
                % e.args)
            return -1.
        else:
            raise e

class muonScaleRes(Module):
    def __init__(self, year, dataset_type, overwritePt=False, maxPt=200., geofit=False):
        self.year = year
        self.dataset_type = dataset_type
        self.overwritePt = overwritePt
        self.maxPt = maxPt
        self.isMC = self.dataset_type.lower() == "mc"
        self.geofit = geofit

        self.era = {"2015": '2016aUL', "2016": '2016bUL', "2017": '2017UL', "2018": '2018UL'}[year]
        correction_file = f"/afs/cern.ch/user/p/pkatris/Run3_Hc/CMSSW_13_3_0/src/PhysicsTools/NanoHc/data/MuonScale/RoccoR{self.era}.txt"

        print(f"Loading muon scale and smearing file from {correction_file}")
        self._roccor = RoccoR(correction_file)
        self.rnd = np.random.default_rng()

    def correct(self, event, mu, variation=None):
        if mu.pt > self.maxPt:
            return mu.pt, mu.pt, mu.pt  # Skip correction for high-pt muons

        roccor = self._roccor
        pt_corr, pt_up, pt_dn = mu.pt, mu.pt, mu.pt
        pt_err = 0

        if not self.isMC:
            pt_corr = mk_safe(roccor.kScaleDT, mu.charge, mu.pt, mu.eta, mu.phi)
            pt_err = mk_safe(roccor.kScaleDTerror, mu.charge, mu.pt, mu.eta, mu.phi)
        else:
            try:
                genparticles = event._genparticles
            except RuntimeError:
                genparticles = Collection(event, "GenPart")
                event._genparticles = genparticles

            if mu.genPartIdx >= 0 and mu.genPartIdx < len(genparticles):
                genMu = genparticles[mu.genPartIdx]
                pt_corr = mk_safe(roccor.k_spread_MC, mu.charge, mu.pt, mu.eta, mu.phi, genMu.pt,0,0)
                pt_err = mk_safe(roccor.k_spread_MC_error, mu.charge, mu.pt, mu.eta, mu.phi, genMu.pt)
            else:
                self.rnd = np.random.default_rng(rndSeed(event, mu))
                u1 = self.rnd.uniform(0.0, 1.0)
                pt_corr = mk_safe(roccor.k_smear_MC, mu.charge, mu.pt, mu.eta, mu.phi, mu.nTrackerLayers, u1,0,0)
                pt_err = mk_safe(roccor.k_smear_MC_error, mu.charge, mu.pt, mu.eta, mu.phi, mu.nTrackerLayers, u1)

        pt_up = pt_corr + pt_err
        pt_dn = pt_corr - pt_err
        # print(pt_corr)
        # print(pt_up)
        # print(pt_dn)

        if self.geofit and abs(mu.dxybs) < 0.01:
            pt_corr = max(0, self._geofit(mu.dxybs * mu.charge, pt_corr, mu.eta, self.year))
            pt_up = max(0, self._geofit(mu.dxybs * mu.charge, pt_up, mu.eta, self.year))
            pt_dn = max(0, self._geofit(mu.dxybs * mu.charge, pt_dn, mu.eta, self.year))

        return max(0, mu.pt * pt_corr), max(0, mu.pt * pt_up), max(0, mu.pt * pt_dn)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if self.overwritePt:
            self.out.branch("Muon_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_uncorrected_pt", "F", lenVar="nMuon")
        else:
            self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_corrected_up_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_corrected_down_pt", "F", lenVar="nMuon")

    def analyze(self, event):
        if event.nMuon == 0:
            return True

        muons = Collection(event, "Muon")
        pt_corr, pt_corr_up, pt_corr_dn = zip(*[self.correct(event, mu) for mu in muons])

        if self.overwritePt:
            self.out.fillBranch("Muon_uncorrected_pt", [mu.pt for mu in muons])
            self.out.fillBranch("Muon_pt", pt_corr)
        else:
            self.out.fillBranch("Muon_corrected_pt", pt_corr)
            self.out.fillBranch("Muon_corrected_up_pt", pt_corr_up)
            self.out.fillBranch("Muon_corrected_down_pt", pt_corr_dn)

        return True
