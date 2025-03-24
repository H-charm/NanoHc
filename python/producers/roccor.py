import os
import math
import numpy as np

class CrystalBall:
    pi = 3.14159
    sqrtPiOver2 = math.sqrt(pi / 2.0)
    sqrt2 = math.sqrt(2.0)

    def __init__(self, s=1, a=10, n=10):
        self.s = s
        self.a = a
        self.n = n
        self.init()

    def init(self):
        fa = abs(self.a)
        ex = math.exp(-fa * fa / 2)
        A = (self.n / fa) ** self.n * ex
        C1 = self.n / fa / (self.n - 1) * ex
        D1 = 2 * CrystalBall.sqrtPiOver2 * math.erf(fa / CrystalBall.sqrt2)

        self.B = self.n / fa - fa
        self.C = (D1 + 2 * C1) / C1
        self.D = (D1 + 2 * C1) / 2
        self.N = 1.0 / self.s / (D1 + 2 * C1)
        self.k = 1.0 / (self.n - 1)

        self.NA = self.N * A
        self.Ns = self.N * self.s
        self.NC = self.Ns * C1
        self.F = 1 - fa * fa / self.n
        self.G = self.s * self.n / fa

    def invcdf(self, u):
        if u < 0.5:
            return self.G * (self.F - (self.NC / u) ** self.k)
        return self.G * (self.F - (self.C - u / self.NC) ** -self.k)

class ResParams:
    def __init__(self):
        self.eta = 0.0
        self.kRes = [1.0, 1.0]  # MC, Data
        self.nTrk = [[], []]  # MC, Data
        self.rsPar = [[], [], []]  # 3 parameters
        self.cb = []

class RocRes:
    MC, Data, Extra = 0, 1, 2
    
    def __init__(self):
        self.NETA = 0
        self.NTRK = 0
        self.NMIN = 0
        self.resol = []
        self.reset()
    
    def reset(self):
        self.NETA = 0
        self.NTRK = 0
        self.NMIN = 0
        self.resol.clear()
    
    def eta_bin(self, eta):
        abseta = abs(eta)
        for i in range(self.NETA - 1):
            if abseta < self.resol[i + 1].eta:
                return i
        return self.NETA - 1
    
    def trk_bin(self, x, h, T=MC):
        for i in range(self.NTRK - 1):
            if x < self.resol[h].nTrk[T][i + 1]:
                return i
        return self.NTRK - 1
    
    def sigma(self, pt, H, F):
        dpt = pt - 45
        rp = self.resol[H]
        return rp.rsPar[0][F] + rp.rsPar[1][F] * dpt + rp.rsPar[2][F] * dpt * dpt
    
    def rndm(self, H, F, w):
        rp = self.resol[H]
        return rp.nTrk[self.MC][F] + (rp.nTrk[self.MC][F + 1] - rp.nTrk[self.MC][F]) * w
    
    def k_spread(self, gpt, rpt, eta, n=None, w=None):
        H = self.eta_bin(abs(eta))
        if n is None:
            k = self.resol[H].kRes
            x = gpt / rpt
            return x / (1.0 + (x - 1.0) * k[self.Data] / k[self.MC])
        else:
            F = max(n - self.NMIN, 0)
            v = self.rndm(H, F, w)
            D = self.trk_bin(v, H, self.Data)
            kold = gpt / rpt
            rp = self.resol[H]
            u = rp.cb[F].cdf((kold - 1.0) / rp.kRes[self.MC] / self.sigma(gpt, H, F))
            knew = 1.0 + rp.kRes[self.Data] * self.sigma(gpt, H, D) * rp.cb[D].invcdf(u)
            return 1.0 if knew < 0 else kold / knew
    
    def k_smear(self, pt, eta, type, v, u, n=None):
        H = self.eta_bin(abs(eta))
        if n is None:
            F = self.trk_bin(v, H)
        else:
            F = n - self.NMIN
            if type == self.Data:
                F = self.trk_bin(self.rndm(H, F, v), H, self.Data)
        rp = self.resol[H]
        x = rp.kRes[type] * self.sigma(pt, H, F) * rp.cb[F].invcdf(u)
        return 1.0 / (1.0 + x)
    
    def k_extra(self, pt, eta, n, u, w=None):
        H = self.eta_bin(abs(eta))
        F = max(n - self.NMIN, 0)
        rp = self.resol[H]
        if w is None:
            d, m = rp.kRes[self.Data], rp.kRes[self.MC]
            x = max(0, np.sqrt(d * d - m * m) * self.sigma(pt, H, F) * rp.cb[F].invcdf(u))
        else:
            v = rp.nTrk[self.MC][F] + (rp.nTrk[self.MC][F + 1] - rp.nTrk[self.MC][F]) * w
            D = self.trk_bin(v, H, self.Data)
            RD = rp.kRes[self.Data] * self.sigma(pt, H, D)
            RM = rp.kRes[self.MC] * self.sigma(pt, H, F)
            x = max(0, np.sqrt(RD * RD - RM * RM) * rp.cb[F].invcdf(u))
        return 1.0 if x <= -1 else 1.0 / (1.0 + x)


class RocOne:
    def __init__(self):
        self.RR = RocRes()
        self.CP = {RocRes.MC: [], RocRes.Data: []}

import numpy as np

class RoccoR:
    MPHI = -np.pi
    
    def __init__(self, filename=None):
        self.NETA = 0
        self.NPHI = 0
        self.etabin = []
        self.nset = 0
        self.nmem = []
        self.RC = []
        if filename:
            self.init(filename)
    
    def reset(self):
        self.NETA = 0
        self.NPHI = 0
        self.etabin.clear()
        self.nset = 0
        self.nmem.clear()
        self.RC.clear()
    
    def init(self, filename):
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        RMIN, RTRK, RETA = 0, 0, 0
        BETA = []
        dKdX = 0
        
        for line in lines:
            tokens = line.split()
            if not tokens:
                continue
            
            if tokens[0] == "VERSION":
                print(f"RoccoR: {tokens[1]}")
                continue
            
            tag = tokens[0]
            if tag == "NSET":
                self.nset = int(tokens[1])
                self.nmem = [0] * self.nset
                self.RC = [[RocOne() for _ in range(self.nmem[i])] for i in range(self.nset)]
            elif tag == "NMEM":
                self.nmem = list(map(int, tokens[1:]))
            elif tag == "RMIN":
                RMIN = int(tokens[1])
            elif tag == "RTRK":
                RTRK = int(tokens[1])
            elif tag == "RETA":
                RETA = int(tokens[1])
                BETA = list(map(float, tokens[2:]))
            elif tag == "CPHI":
                self.NPHI = int(tokens[1])
                self.DPHI = 2 * np.pi / self.NPHI
            elif tag == "CETA":
                self.NETA = int(tokens[1])
                self.etabin = list(map(float, tokens[2:]))
    
    def eta_bin(self, x):
        for i in range(self.NETA - 1):
            if x < self.etabin[i + 1]:
                return i
        return self.NETA - 1
    
    def phi_bin(self, x):
        ibin = int((x - self.MPHI) / self.DPHI)
        return max(0, min(ibin, self.NPHI - 1))
    
    def k_scale_DT(self, Q, pt, eta, phi, s, m):
        H, F = self.eta_bin(eta), self.phi_bin(phi)
        return self.RC[s][m].CP[RocRes.Data][H][F].k(Q, pt)
    
    def k_scale_MC(self, Q, pt, eta, phi, s, m):
        H, F = self.eta_bin(eta), self.phi_bin(phi)
        return self.RC[s][m].CP[RocRes.MC][H][F].k(Q, pt)
    
    def k_spread_MC(self, Q, pt, eta, phi, gt, s, m):
        rc = self.RC[s][m]
        H, F = self.eta_bin(eta), self.phi_bin(phi)
        k = rc.CP[RocRes.MC][H][F].k(Q, pt)
        return k * rc.RR.k_spread(gt, k * pt, eta)
    
    def k_smear_MC(self, Q, pt, eta, phi, n, u, s, m):
        rc = self.RC[s][m]
        H, F = self.eta_bin(eta), self.phi_bin(phi)
        k = rc.CP[RocRes.MC][H][F].k(Q, pt)
        return k * rc.RR.k_extra(k * pt, eta, n, u)
    
    def k_gen_smear(self, pt, eta, v, u, TT, s, m):
        if not self.RC:
            return 1.0
        return self.RC[s][m].RR.kSmear(pt, eta, TT, v, u)
    
    def error(self, func):
        sum_sq = sum((func(s, i) - func(0, 0)) ** 2 / self.nmem[s] for s in range(self.nset) for i in range(self.nmem[s]))
        return np.sqrt(sum_sq)
    
    def k_scale_DT_error(self, Q, pt, eta, phi):
        return self.error(lambda s, m: self.k_scale_DT(Q, pt, eta, phi, s, m))
    
    def k_spread_MC_error(self, Q, pt, eta, phi, gt):
        return self.error(lambda s, m: self.k_spread_MC(Q, pt, eta, phi, gt, s, m))
    
    def k_smear_MC_error(self, Q, pt, eta, phi, n, u):
        return self.error(lambda s, m: self.k_smear_MC(Q, pt, eta, phi, n, u, s, m))
