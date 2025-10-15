# zz_kfactors_apply_by_sample.py
import math
import json
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

# -----------------------------
# Utilities
# -----------------------------
def _finite(x, default=1.0):
    try:
        xf = float(x)
        return xf if math.isfinite(xf) else float(default)
    except Exception:
        return float(default)

def _evalSpline(sp, x):
    if not sp:
        return 1.0
    xmin = sp.GetXmin(); xmax = sp.GetXmax()
    if x < xmin: x = xmin
    if x > xmax: x = xmax
    return sp.Eval(x)

def _derive_from_best(event):
    try:
        hcands = getattr(event, "Hcandidates", [])
        if not hcands:
            return 0.0, 1
        best = hcands[0]
        mzz = float(best.mass)
        ids = {
            abs(best.Z1.lep1.pdgId), abs(best.Z1.lep2.pdgId),
            abs(best.Z2.lep1.pdgId), abs(best.Z2.lep2.pdgId)
        }
        flavor = 1 if (ids == {11} or ids == {13}) else 2
        return mzz, flavor
    except Exception:
        return 0.0, 1

def _pt_dphi_from_best(event):
    try:
        hcands = getattr(event, "Hcandidates", [])
        if not hcands:
            return 0.0, 0.0
        best = hcands[0]
        Z1 = getattr(best, "Z1", None)
        Z2 = getattr(best, "Z2", None)

        def _phi(o):
            for attr in ("phi", "Phi"):
                if hasattr(o, attr):
                    v = getattr(o, attr)
                    return float(v() if callable(v) else v)
            if hasattr(o, "p4"):
                p4 = o.p4() if callable(o.p4) else o.p4
                if hasattr(p4, "Phi"):
                    return float(p4.Phi())
            return 0.0

        def _pt(o):
            for attr in ("pt", "Pt"):
                if hasattr(o, attr):
                    v = getattr(o, attr)
                    return float(v() if callable(v) else v)
            if hasattr(o, "p4"):
                p4 = o.p4() if callable(o.p4) else o.p4
                if hasattr(p4, "Pt"):
                    return float(p4.Pt())
            return 0.0

        phi1 = _phi(Z1) if Z1 is not None else 0.0
        phi2 = _phi(Z2) if Z2 is not None else 0.0
        dphi = abs(ROOT.TVector2.Phi_mpi_pi(phi1 - phi2))

        pt1 = _pt(Z1) if Z1 is not None else 0.0
        pt2 = _pt(Z2) if Z2 is not None else 0.0
        ptZZ = math.sqrt(
            (pt1*math.cos(phi1) + pt2*math.cos(phi2))**2 +
            (pt1*math.sin(phi1) + pt2*math.sin(phi2))**2
        )
        return float(ptZZ), float(dphi)
    except Exception:
        return 0.0, 0.0


# -----------------------------
# GG → ZZ k-factor
# -----------------------------
class GGZZKFactorProducer(Module):
    def __init__(self,
                 year,
                 dataset_type,
                 gg_nnlo_file=None,
                 gg_nlo_file=None,
                 gg_mode='NNLO_NLO',
                 sample_name=None,             # current sample name
                 apply_for_samples=None,       # list of names to apply
                 debug=False):
        assert gg_mode in ('NNLO_NLO', 'NNLO_LO', 'NLO_LO')
        self.year = year
        self.dataset_type = dataset_type
        self.gg_nnlo_file = gg_nnlo_file
        self.gg_nlo_file  = gg_nlo_file
        self.gg_mode = gg_mode
        self.sample_name = str(sample_name) if sample_name else None
        self.apply_for_samples = set(apply_for_samples or [])
        self.debug = bool(debug)

        self._sp_ggNNLO = None
        self._sp_ggNLO  = None
        self._warned_once = False
        self._apply = False
        self.isMC = False
        self.out = None
        self.branch_name = 'KFactor_QCD_ggZZ'

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = (self.dataset_type == "mc")
        if not self.isMC:
            return

        # check if this sample should get the weight
        self._apply = (self.sample_name in self.apply_for_samples)

        self.out = wrappedOutputTree
        self.out.branch(self.branch_name, "F")

        if not self._apply:
            return  # skip spline loading

        if self.gg_nnlo_file:
            fNNLO = ROOT.TFile.Open(self.gg_nnlo_file)
            sp = fNNLO.Get("sp_kfactor_Nominal")
            self._sp_ggNNLO = sp
            fNNLO.Close()
        if self.gg_nlo_file:
            fNLO = ROOT.TFile.Open(self.gg_nlo_file)
            sp = fNLO.Get("sp_kfactor_Nominal")
            self._sp_ggNLO = sp
            fNLO.Close()

    def analyze(self, event):
        if not self.isMC:
            return True

        k_qcd_gg = 1.0
        if self._apply:
            mZZ, _ = _derive_from_best(event)
            if mZZ > 0.0:
                try:
                    if self.gg_mode == 'NNLO_NLO' and self._sp_ggNNLO and self._sp_ggNLO:
                        num = _evalSpline(self._sp_ggNNLO, mZZ)
                        den = _evalSpline(self._sp_ggNLO , mZZ)
                        k_qcd_gg = float(num / den) if den else 1.0
                    elif self.gg_mode == 'NNLO_LO' and self._sp_ggNNLO:
                        k_qcd_gg = float(_evalSpline(self._sp_ggNNLO, mZZ))
                    elif self.gg_mode == 'NLO_LO' and self._sp_ggNLO:
                        k_qcd_gg = float(_evalSpline(self._sp_ggNLO, mZZ))
                except Exception as e:
                    if not self._warned_once:
                        print(f"[GGZZKFactorProducer] WARNING: eval failed ({e}); using 1.0")
                        self._warned_once = True

        self.out.fillBranch(self.branch_name, _finite(k_qcd_gg, 1.0))
        return True


# -----------------------------
# QQ → ZZ k-factor (JSON only)
# -----------------------------
class QQZZKFactorProducer(Module):
    def __init__(self,
                 year,
                 dataset_type,
                 json_table_path,
                 write_pt_branch=False,
                 write_dphi_branch=False,
                 sample_name=None,             # current sample name
                 apply_for_samples=None,       # list of names to apply
                 debug=False,
                 print_every=1000):
        self.year = year
        self.dataset_type = dataset_type
        self.json_table_path = json_table_path
        self.write_pt_branch = bool(write_pt_branch)
        self.write_dphi_branch = bool(write_dphi_branch)
        self.sample_name = str(sample_name) if sample_name else None
        self.apply_for_samples = set(apply_for_samples or [])
        self.debug = bool(debug)
        self.print_every = int(print_every)

        self._tbl = None
        self._apply = False
        self._iev = 0
        self.isMC = False
        self.out = None
        self._warned_once = False

    def _load_json(self, path):
        with open(path) as f:
            return json.load(f)

    def _kfactor_M_json(self, mzz, finalState, order):
        rows = self._tbl["mass_xsec"][str(finalState)]
        edges = [r[0] for r in rows]
        i = 0
        for ii in range(1, len(edges)):
            if mzz < edges[ii]:
                i = ii - 1
                break
        else:
            i = len(rows) - 1
        _, lo, nlo, nnlo = rows[i]
        num = nlo if order == 1 else nnlo
        return (num / lo) if lo > 0 else 1.0

    def _kfactor_Pt_json(self, ptzz, finalState):
        for lo, hi, val in self._tbl["pt_bins"][str(finalState)]:
            if (ptzz > lo) and (ptzz <= hi):
                return float(val)
        if ptzz > self._tbl["pt_bins"][str(finalState)][-1][1]:
            return float(self._tbl["pt_overflow"][str(finalState)])
        return 1.0

    def _kfactor_dPhi_json(self, abs_dphi, finalState):
        for lo, hi, val in self._tbl["dphi_bins"][str(finalState)]:
            if (abs_dphi > lo) and (abs_dphi <= hi):
                return float(val)
        if abs_dphi > 2.9:
            return float(self._tbl["dphi_overflow"][str(finalState)])
        return 1.0

    def beginJob(self):
        try:
            self._tbl = self._load_json(self.json_table_path)
            print(f"[QQZZKFactorProducer] Loaded JSON {self.json_table_path}")
        except Exception as e:
            self._tbl = None
            print(f"[QQZZKFactorProducer] WARNING: JSON load failed ({e}); set to 1.0")

        # check if this sample should get the weight
        self._apply = (self.sample_name in self.apply_for_samples)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = (self.dataset_type == "mc")
        if not self.isMC:
            return
        self.out = wrappedOutputTree

        self.out.branch("KFactor_QCD_qqZZ_M", "F")
        if self.write_pt_branch:
            self.out.branch("KFactor_QCD_qqZZ_Pt", "F")
        if self.write_dphi_branch:
            self.out.branch("KFactor_QCD_qqZZ_dPhi", "F")

    def analyze(self, event):
        if not self.isMC:
            return True
        self._iev += 1

        k_M = 1.0
        k_Pt = 1.0
        k_dPhi = 1.0

        if self._apply and (self._tbl is not None):
            mZZ, flavor = _derive_from_best(event)
            ptZZ, dphiZZ = _pt_dphi_from_best(event)

            if mZZ > 0.0 and (flavor in (1, 2)):
                try:
                    k_nlo_over_lo  = self._kfactor_M_json(mZZ, flavor, 1)
                    k_nnlo_over_lo = self._kfactor_M_json(mZZ, flavor, 2)
                    k_M = float(k_nnlo_over_lo / k_nlo_over_lo) if k_nlo_over_lo else 1.0
                except Exception as e:
                    if not self._warned_once:
                        print(f"[QQZZKFactorProducer] WARNING: mZZ eval failed ({e}); using 1.0")
                        self._warned_once = True

            if self.write_pt_branch and (ptZZ > 0.0):
                try:
                    k_Pt = float(self._kfactor_Pt_json(ptZZ, flavor))
                except Exception:
                    pass

            if self.write_dphi_branch and (dphiZZ > 0.0):
                try:
                    k_dPhi = float(self._kfactor_dPhi_json(dphiZZ, flavor))
                except Exception:
                    pass

        self.out.fillBranch("KFactor_QCD_qqZZ_M", _finite(k_M, 1.0))
        if self.write_pt_branch:
            self.out.fillBranch("KFactor_QCD_qqZZ_Pt", _finite(k_Pt, 1.0))
        if self.write_dphi_branch:
            self.out.fillBranch("KFactor_QCD_qqZZ_dPhi", _finite(k_dPhi, 1.0))

        return True


