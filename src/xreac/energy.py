"""Differentiable Zn/O ReaxFF energies in kcal/mol and Angstrom.

Equations and conventions adapted from the LAMMPS/PuReMD REAXFF sources,
stable_22Jul2025_update4. Copyright (2010) Purdue University; LAMMPS
contributors. Distributed under GPL-2.0-or-later; see NOTICE and LICENSE.
"""
from itertools import combinations

import autograd.numpy as np
import numpy as onp
from autograd.tracer import getval

C_ELE = 332.06371
QEQ_COULOMB = 14.4
SELF_CONVERSION = 23.02
THB_CUT = .001
THB_PRODUCT_CUT = .00001
BOND_CUT = 5.0
COMPONENTS = ("bond", "atom", "lone_pair", "molecule", "angle", "penalty",
              "angle_conjugation", "hydrogen_bond", "torsion", "conjugation",
              "vdw", "coulomb", "electric_field", "qeq")


def positive_power(x, p):
    # Avoid undefined derivatives of masked zero raised to fractional powers.
    return np.where(x > 0, np.where(x > 0, x, 1.0) ** p, 0.0)


class EnergyModel:
    def __init__(self, ff, symbols):
        self.ff = ff
        self.symbols = tuple(symbols)
        self.n = len(symbols)
        self.g = ff.general
        self.a = {k: onp.array([ff.atoms[s][k] for s in symbols])
                  for k in ff.atoms[symbols[0]]}
        self.p = {k: onp.array([[ff.pairs[s, t][k] for t in symbols] for s in symbols])
                  for k in ff.pairs[symbols[0], symbols[0]]}

    def geometry(self, x):
        delta = x[:, None, :] - x[None, :, :]
        r = np.sqrt(np.sum(delta * delta, axis=-1) + np.eye(self.n))
        return delta, r

    def electrostatics(self, r):
        u = np.minimum(r / self.g[12], 1.0)
        taper = 1 + u**4 * (-35 + u * (84 + u * (-70 + 20*u)))
        taper = np.where(r < self.g[12], taper, 0.0) * (1-np.eye(self.n))
        shield = taper / (r**3 + self.p["gamma"]) ** (1/3)
        h = QEQ_COULOMB * shield + np.diag(self.a["eta"])
        # A Lagrange multiplier enforces total charge = 0 exactly.
        ones = np.ones((self.n, 1))
        kkt = np.concatenate((np.concatenate((h, ones), axis=1),
                              np.concatenate((ones.T, np.zeros((1, 1))), axis=1)), axis=0)
        rhs = np.concatenate((-self.a["chi"], np.zeros(1)))
        q = np.linalg.solve(kkt, rhs)[:self.n]
        return q, taper, shield

    def bond_orders(self, r):
        p, a, g = self.p, self.a, self.g
        cutoff = .01 * g[29]
        sigma = (1+cutoff) * np.exp(p["p_bo1"] * (r/p["r_s"])**p["p_bo2"])
        pi_ok = (a["r_pi"][:, None] > 0) & (a["r_pi"][None, :] > 0)
        pp_ok = (a["r_pi_pi"][:, None] > 0) & (a["r_pi_pi"][None, :] > 0)
        pi = np.where(pi_ok, np.exp(p["p_bo3"] * (r / np.where(pi_ok, p["r_p"], 1))**p["p_bo4"]), 0)
        pp = np.where(pp_ok, np.exp(p["p_bo5"] * (r / np.where(pp_ok, p["r_pp"], 1))**p["p_bo6"]), 0)
        raw = sigma + pi + pp
        mask = (getval(raw) >= cutoff) & (getval(r) <= BOND_CUT) & ~onp.eye(self.n, dtype=bool)
        bo = np.where(mask, raw-cutoff, 0)
        pi, pp = np.where(mask, pi, 0), np.where(mask, pp, 0)
        total = np.sum(bo, axis=1)
        d = total-a["valency"]
        db = total-a["valency_boc"]
        f2 = np.exp(-g[0]*d[:, None]) + np.exp(-g[0]*d[None, :])
        f3 = -np.log(.5*(np.exp(-g[1]*d[:, None])+np.exp(-g[1]*d[None, :]))) / g[1]
        vi, vj = a["valency"][:, None], a["valency"][None, :]
        f1 = np.where(p["ovc"] >= .001, .5*((vi+f2)/(vi+f2+f3)+(vj+f2)/(vj+f2+f3)), 1)
        f4 = 1/(1+np.exp(-(p["p_boc4"]*bo**2-db[:, None])*p["p_boc3"]+p["p_boc5"]))
        f5 = 1/(1+np.exp(-(p["p_boc4"]*bo**2-db[None, :])*p["p_boc3"]+p["p_boc5"]))
        corr = f1 * np.where(p["v13cor"] >= .001, f4*f5, 1)
        bo, pi, pp = bo*corr, pi*corr*f1, pp*corr*f1
        sigma = bo-pi-pp
        return tuple(np.where(v >= 1e-10, v, 0.0) for v in (bo, sigma, pi, pp))

    def components(self, x, fixed_charges=None):
        p, a, g = self.p, self.a, self.g
        delta, r = self.geometry(x)
        q, taper, shield = self.electrostatics(r)
        if fixed_charges is not None:
            q = fixed_charges
        bo, sigma, pi, pp = self.bond_orders(r)
        total = np.sum(bo, axis=1)
        d, db = total-a["valency"], total-a["valency_boc"]
        de = total-a["valency_e"]
        half = onp.trunc(getval(de)/2)
        vlpex = de-2*half
        nlp = np.exp(-g[15]*(2+vlpex)**2)-half
        dlp = .5*(a["valency_e"]-a["valency"])-nlp
        dlpt = np.where(a["mass"] > 21, 0, dlp)
        e = {k: 0.0 for k in COMPONENTS}
        e["bond"] = -.5*np.sum(p["De_s"]*sigma*np.exp(p["p_be1"]*(1-positive_power(sigma, p["p_be2"])))
                              +p["De_p"]*pi+p["De_pp"]*pp)
        e["lone_pair"] = np.sum(a["p_lp2"]*dlp/(1+np.exp(-75*dlp)))
        dfvl = (a["mass"] <= 21).astype(float)
        ov1 = np.sum(p["p_ovun1"]*p["De_s"]*bo, axis=1)
        ov2 = np.sum((d[None, :]-dfvl[:, None]*dlpt[None, :])*(pi+pp), axis=1)
        dc = d-dfvl*dlpt/(1+g[32]*np.exp(g[31]*ov2))
        over = ov1*dc/(dc+a["valency"]+1e-8)/(1+np.exp(a["p_ovun2"]*dc))
        under = -a["p_ovun5"]*(1-np.exp(g[6]*dc))/(1+np.exp(-a["p_ovun2"]*dc))/(1+g[8]*np.exp(g[9]*ov2))
        e["atom"] = np.sum(over+under)
        f13 = (r**g[28]+(1/p["gamma_w"])**g[28])**(1/g[28])
        ev = np.exp(.5*p["alpha"]*(1-f13/p["r_vdW"]))
        e["vdw"] = .5*np.sum(taper*p["D"]*(ev**2-2*ev))
        e["coulomb"] = .5*C_ELE*np.sum(q[:, None]*q[None, :]*shield)
        e["qeq"] = SELF_CONVERSION*np.sum(a["chi"]*q+.5*a["eta"]*q*q)
        neighbors = [onp.flatnonzero(row > THB_CUT) for row in getval(bo)]
        self.angles(e, delta, r, bo, pi, pp, total, d, db, nlp, vlpex, neighbors)
        self.torsions(e, delta, r, bo, pi, db, neighbors)
        return np.stack([e[k] for k in COMPONENTS]), q

    def angles(self, e, delta, r, bo, pi, pp, total, d, db, nlp, vlpex, neighbors):
        a, g = self.a, self.g
        triples, params = [], []
        bv = getval(bo)
        for j, ns in enumerate(neighbors):
            for i, k in combinations(ns, 2):
                if bv[i, j]*bv[j, k] <= THB_PRODUCT_CUT:
                    continue
                for prm in self.ff.angles.get((self.symbols[i], self.symbols[j], self.symbols[k]), []):
                    if abs(prm[1]) > .001:
                        triples.append((i, j, k))
                        params.append(prm)
        if not triples:
            return
        i, j, k = onp.array(triples).T
        theta00, v1, v2, coa1, v7, pen1, v4 = onp.array(params).T
        prod = np.exp(-np.sum(bo**8, axis=1))
        sbo = np.sum(pi+pp, axis=1)+(1-prod)*(-db-g[33]*np.where(getval(vlpex) >= 0, 0, nlp))
        sbo2 = np.where(sbo <= 0, 0, np.where(sbo <= 1, positive_power(sbo, g[16]),
                        np.where(sbo < 2, 2-positive_power(2-sbo, g[16]), 2)))
        theta0 = (180-theta00*(1-np.exp(-g[17]*(2-sbo2[j]))))*(3.14159265/180)
        cos = np.sum(delta[i, j]*delta[k, j], axis=1)/(r[i, j]*r[k, j])
        # Exact linear angles are nondifferentiable; keep arccos finite there.
        theta = np.arccos(np.clip(cos, -1+1e-14, 1-1e-14))
        bij, bjk = bo[i, j]-THB_CUT, bo[j, k]-THB_CUT
        f7 = (1-np.exp(-a["p_val3"][j]*bij**v4))*(1-np.exp(-a["p_val3"][j]*bjk**v4))
        ex6, ex7 = np.exp(g[14]*db[j]), np.exp(-v7*db[j])
        f8 = a["p_val5"][j]-(a["p_val5"][j]-1)*(2+ex6)/(1+ex6+ex7)
        exangle = np.exp(-v2*(theta0-theta)**2)
        e["angle"] = np.sum(f7*f8*v1*np.where(v1 >= 0, 1-exangle, -exangle))
        ex3, ex4 = np.exp(-g[20]*d[j]), np.exp(g[21]*d[j])
        e["penalty"] = np.sum(pen1*(2+ex3)/(1+ex3+ex4)*np.exp(-g[19]*((bij-2)**2+(bjk-2)**2)))
        e["angle_conjugation"] = np.sum(coa1/(1+np.exp(g[2]*(total[j]-a["valency_val"][j])))
            *np.exp(-g[38]*((total[i]-bij)**2+(total[k]-bjk)**2))
            *np.exp(-g[30]*((bij-1.5)**2+(bjk-1.5)**2)))

    def torsions(self, e, delta, r, bo, pi, db, neighbors):
        quads, params = [], []
        bv = getval(bo)
        for j, ns in enumerate(neighbors):
            for k in ns:
                if j >= k:
                    continue
                for i in ns:
                    if i == k:
                        continue
                    for l in neighbors[k]:
                        if l == j or l == i or bv[i, j]*bv[j, k]*bv[k, l] <= THB_CUT:
                            continue
                        prm = self.ff.torsions.get(tuple(self.symbols[t] for t in (i, j, k, l)))
                        if prm is not None:
                            quads.append((i, j, k, l))
                            params.append(prm)
        if not quads:
            return
        i, j, k, l = onp.array(quads).T
        v1, v2, v3, tor1, cot1 = onp.array(params).T
        # LAMMPS vectors point j->i, j->k, k->l.
        u, v, w = delta[i, j], delta[k, j], delta[l, k]
        cos1 = np.sum(u*v, axis=1)/(r[i, j]*r[j, k])
        cos2 = -np.sum(v*w, axis=1)/(r[j, k]*r[k, l])
        sinprod = np.sqrt(np.maximum(1-cos1*cos1, 1e-20))*np.sqrt(np.maximum(1-cos2*cos2, 1e-20))
        normcos = r[j, k]**2*np.sum(u*w, axis=1)-np.sum(u*v, axis=1)*np.sum(v*w, axis=1)
        normsin = -r[j, k]*np.sum(u*np.cross(v, w), axis=1)
        cw = normcos/np.sqrt(np.maximum(normcos**2+normsin**2, 1e-30))
        b1, b2, b3 = bo[i, j]-THB_CUT, bo[j, k]-THB_CUT, bo[k, l]-THB_CUT
        ds = db[j]+db[k]
        ex3, ex4 = np.exp(-self.g[24]*ds), np.exp(self.g[25]*ds)
        f11 = (2+ex3)/(1+ex3+ex4)
        ex1 = np.exp(tor1*(2-pi[j, k]-f11)**2)
        fn10 = (1-np.exp(-self.g[23]*b1))*(1-np.exp(-self.g[23]*b2))*(1-np.exp(-self.g[23]*b3))
        cv = .5*(v1*(1+cw)+v2*ex1*(2-2*cw*cw)+v3*(1+4*cw**3-3*cw))
        e["torsion"] = np.sum(fn10*sinprod*cv)
        fn12 = np.exp(-self.g[27]*((b1-1.5)**2+(b2-1.5)**2+(b3-1.5)**2))
        e["conjugation"] = np.sum(cot1*fn12*(1+(cw*cw-1)*sinprod))
