import numpy as np
import os
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

era_dict = {"2022": '2022_Summer22', "2022EE": '2022_Summer22EE', "2023": '2023_Summer23', "2023BPix": '2023_Summer23BPix'}

# year: L1Key, L2Key, L3Key, L2L3Key, JERKey. JERsfKey
keyMC_dict = {
    "2022": [
        "Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi",
        "Summer22_22Sep2023_V2_MC_L2Relative_AK4PFPuppi",
        "Summer22_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi",
        "Summer22_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi",
        "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
        "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
    ],
    "2022EE": [
        "Summer22EE_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_V2_MC_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi",
        "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
        "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
    ],
    "2023": [
        "Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_MC_L2L3Residual_AK4PFPuppi",
        "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi",
        "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"
    ],
    "2023BPix": [
        "Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi",
        "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi",
        "Summer23BPixPrompt23_V3_MC_L3Absolute_AK4PFPuppi",
        "Summer23BPixPrompt23_V3_MC_L2L3Residual_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
    ]
}

# year: L1Key, L2Key, L3Key, L2L3Key
keyData_dict={
    # 2022
    "Run2022C-22Sep2023-v1": [
        "Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2022D-22Sep2023-v1": [
        "Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22_22Sep2023_RunCD_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    # 2022EE
    #eraE
    "Run2022E-22Sep2023-v1": [
        "Summer22EE_22Sep2023_RunE_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunE_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunE_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunE_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    #eraF
    "Run2022F-22Sep2023-v1": [
        "Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2022F-22Sep2023-v2":[
        "Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunF_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    #eraG
    "Run2022G-22Sep2023-v1": [
        "Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2022G-22Sep2023-v2": [
        "Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer22EE_22Sep2023_RunG_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    # 2023
    "Run2023C-22Sep2023_v1-v1": [
        "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2023C-22Sep2023_v2-v1": [
        "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2023C-22Sep2023_v3-v1": [
        "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2023C-22Sep2023_v4-v1": [
        "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2023C-22Sep2023_v4-v2": [
        "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi",
        "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
    ],
    #2023BPix
    "Run2023D-22Sep2023_v1-v1": [
        "Summer23BPixPrompt23_RunD_V3_DATA_L1FastJet_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L2Relative_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L3Absolute_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L2L3Residual_AK4PFPuppi"
    ],
    "Run2023D-22Sep2023_v2-v1": [
        "Summer23BPixPrompt23_RunD_V3_DATA_L1FastJet_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L2Relative_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L3Absolute_AK4PFPuppi",
        "Summer23BPixPrompt23_RunD_V3_DATA_L2L3Residual_AK4PFPuppi"
    ]
}
class JetJERCProducer(Module):
    def __init__(self, year, era_data, dataset_type, overwritePt=True, usePhiDependentJEC=False, **kwargs):
        """Correct jets following recommendations of JME POG.
        Parameters:
            json_JERC: full path of json file with JERC corrections
            json_JERsmear: full path of json file with smearing terms for JER
            L1Key: key for L1 corrections
            L2Key: key for L2 corrections
            L3Key: key for L3 corrections
            L2L3Key: key for residual corrections
            smearKey: key for smearing formula (None for Data)
            JERKey: key for JER (None for Data)
            JERsfKey: key for JER scale factor (None for Data)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """

        self.year = year
        self.era_data = era_data
        self.dataset_type = dataset_type
        self.era = era_dict[self.year]
        self.isMC = True if self.dataset_type == "mc" else False

        correction_file_JERC = f'/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{self.era}/jet_jerc.json.gz'
        correction_file_JERsmear = '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz'

        self.overwritePt = overwritePt
        self.usePhiDependentJEC = True if year=="2023BPix" else False

        self.evaluator_JERC = correctionlib.CorrectionSet.from_file(correction_file_JERC)
        self.evaluator_jer = correctionlib.CorrectionSet.from_file(correction_file_JERsmear)

        if self.isMC:
            ## Starting from Run 3, the use of PUPPI jets eliminates the need for L1 corrections. To maintain compatibility with Run 2 scripts, a dummy file is provided.
            self.evaluator_L1 = self.evaluator_JERC[keyMC_dict[self.year][0]]
            ## For the next step, there is a mismatch between the terminology in the twiki and the tags in correctionlib
            ## MC-truth = L2 + L3
            self.evaluator_L2 = self.evaluator_JERC[keyMC_dict[self.year][1]]
            self.evaluator_L3 = self.evaluator_JERC[keyMC_dict[self.year][2]]
            ## L2L3residuals should be applied only to data
            ## A dummy file with all entries equal to 1 is provided for MC to have a common script for both data and MC
            self.evaluator_L2L3 = self.evaluator_JERC[keyMC_dict[self.year][3]]
            self.evaluator_JERsmear = self.evaluator_jer['JERSmear']
            self.evaluator_JER = self.evaluator_JERC[keyMC_dict[self.year][4]]
            self.evaluator_JERsf = self.evaluator_JERC[keyMC_dict[self.year][5]]
        else:
            self.evaluator_L1 = self.evaluator_JERC[keyData_dict[self.era_data][0]]
            self.evaluator_L2 = self.evaluator_JERC[keyData_dict[self.era_data][1]]
            self.evaluator_L3 = self.evaluator_JERC[keyData_dict[self.era_data][2]]
            self.evaluator_L2L3 = self.evaluator_JERC[keyData_dict[self.era_data][3]]
            

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Jet_pt", "F", lenVar="nJet")
            self.out.branch("Jet_mass", "F", lenVar="nJet")
            self.out.branch("Jet_uncorrected_pt", "F", lenVar="nJet")
            self.out.branch("Jet_uncorrected_mass", "F", lenVar="nJet")
        else:
            self.out.branch("Jet_corrected_pt", "F", lenVar="nJet")
            self.out.branch("Jet_corrected_mass", "F", lenVar="nJet")
        if self.isMC:
            self.out.branch("Jet_scaleUp_pt", "F", lenVar="nJet")
            self.out.branch("Jet_scaleDn_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearUp_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearDn_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearUp_mass", "F", lenVar="nJet")
            self.out.branch("Jet_smearDn_mass", "F", lenVar="nJet")

    def fixPhi(self, phi):
        if phi > np.pi:
            phi -= 2*np.pi
        elif phi < -np.pi:
            phi += 2*np.pi
        return phi

    def analyze(self, event):
        jets = Collection(event, "Jet")
        if self.isMC: ## genJet info is necessary for JER
            gen_jets = Collection(event, "GenJet")
            gen_jets_pt = np.array([gen_jet.pt for gen_jet in gen_jets])
            gen_jets_eta = np.array([gen_jet.eta for gen_jet in gen_jets])
            gen_jets_phi = np.array([gen_jet.phi for gen_jet in gen_jets])

        pt_corr = []
        pt_uncorr = []
        mass_corr = []
        mass_uncorr = []
        pt_smear_up = []
        pt_smear_dn = []
        mass_smear_up = []
        mass_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []

        for jet in jets:
            #### JEC ####
            ## To be applied to both data and MC
            ## Jet in NanoAOD are already corrected
            ## The correction should be removed and the latest one available should be applied
            pt_raw = jet.pt * (1 - jet.rawFactor)
            mass_raw = jet.mass * (1 - jet.rawFactor)
            ## The three steps of JEC corrections are provided separately
            pt_L1 = pt_raw * self.evaluator_L1.evaluate(jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll)

            if self.usePhiDependentJEC:
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, jet.phi, pt_L1)
            else:
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, pt_L1)

            pt_L3 = pt_L2 * self.evaluator_L3.evaluate(jet.eta, pt_L2)
            if self.year == "2023":
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(float(event.run), jet.eta, pt_L3)
            else:
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(jet.eta, pt_L3)
            JEC = pt_JEC / pt_raw
            mass_JEC = mass_raw * JEC

            if self.isMC:
                #### JER ####
                ## Hybrid method is implemented [https://cms-jerc.web.cern.ch/JER/#smearing-procedures]
                JER = self.evaluator_JER.evaluate(jet.eta, pt_JEC, event.Rho_fixedGridRhoFastjetAll)
                ## GenMatching with genJet
                delta_eta = jet.eta - gen_jets_eta
                fixPhi = np.vectorize(self.fixPhi, otypes=[float])
                delta_phi = fixPhi(jet.phi - gen_jets_phi)
                pt_gen = np.where((np.abs(pt_JEC - gen_jets_pt) < 3 * pt_JEC * JER) & (np.sqrt(delta_eta**2 + delta_phi**2)<0.2), gen_jets_pt, -1.0)
                pt_gen = pt_gen[pt_gen > 0][0] if np.any(pt_gen > 0) else -1. ## If no gen-matching, simply -1
                JERsf = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "nom")
                JERsf_up = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "up")
                JERsf_dn = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "down")
                JERsmear = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf)
                JERsmear_up = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_up)
                JERsmear_dn = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_dn)
                pt_JEC_JER = pt_JEC * JERsmear
                pt_JEC_JER_up = pt_JEC * JERsmear_up
                pt_JEC_JER_dn = pt_JEC * JERsmear_dn
                mass_JEC_JER = mass_JEC * JERsmear
                mass_JEC_JER_up = mass_JEC * JERsmear_up
                mass_JEC_JER_dn = mass_JEC * JERsmear_dn

                pt_corr.append(pt_JEC_JER)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC_JER)
                mass_uncorr.append(mass_raw)
                pt_smear_up.append(pt_JEC_JER_up)
                pt_smear_dn.append(pt_JEC_JER_dn)
                mass_smear_up.append(mass_JEC_JER_up)
                mass_smear_dn.append(mass_JEC_JER_dn)
                pt_scale_up.append(-1)
                pt_scale_dn.append(-1)
            else:
                ## Data
                ## No JER for Data
                pt_corr.append(pt_JEC)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC)
                mass_uncorr.append(mass_raw)

        if self.overwritePt :
            self.out.fillBranch("Jet_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Jet_pt", pt_corr)
            self.out.fillBranch("Jet_uncorrected_mass", mass_uncorr)
            self.out.fillBranch("Jet_mass", mass_corr)
        else :
            self.out.fillBranch("Jet_corrected_pt", pt_corr)
            self.out.fillBranch("Jet_corrected_mass", mass_corr)

        if self.isMC:
            self.out.fillBranch("Jet_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Jet_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Jet_smearUp_mass", mass_smear_up)
            self.out.fillBranch("Jet_smearDn_mass", mass_smear_dn)
            self.out.fillBranch("Jet_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Jet_scaleDn_pt", pt_scale_dn)

        return True